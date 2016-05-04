#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <sys/time.h>
#include <time.h>

#include "cachesim.h"

//#define debug
//#define print_file

//Function prototypes
struct cache_struct** allocate_cache(uint16_t CACHE_DEPTH, uint16_t N_WAY);
void initialize_cache(void);
bool check_hit_cache(uint64_t compare_tag, uint16_t index, uint16_t *way);
bool check_cache_full(uint16_t index, uint16_t *empty_way);
bool check_cache_empty(uint16_t index, uint16_t *filled_way);
uint16_t LRU_way_cache(uint16_t index);

struct VC_struct* allocate_VC(uint8_t VC_DEPTH);
void initialize_VC(void);
bool check_hit_VC(uint64_t compare_tag, uint16_t cache_index, uint8_t *VC_row);
bool check_VC_full(uint8_t *empty_row);
bool check_VC_empty(uint8_t *filled_row);
uint8_t first_in_row_VC(void);

void prefetch(uint64_t fetch_address, struct cache_stats_t* p_stats);

void Get_tag_index_offset(uint64_t address, uint64_t *tag, uint16_t *index, uint16_t *offset);

// Global variables
uint8_t C;
uint8_t B;
uint8_t S;
uint16_t CACHE_DEPTH; // Cache depth
uint16_t N_WAY; // No. of ways
uint8_t VC_DEPTH;
uint8_t PREFETCH_DEGREE;

uint64_t Last_Miss_Block_Addr;
__int128_t pending_stride; // Pending stride can be negative.

#ifdef print_file
  FILE *fw;
#endif


// If way-0 of cache is empty then entire cache is empty
// If row-0 of VC is empty then entire VC is empty

// This structure template represents 1-cache-line.
struct cache_struct
{
  uint64_t tag;
  bool valid;
  bool dirty;
  long last_acc_time; // Time stamp indicating the last accessed time
  bool prefetch_not_acc; // 'false' if miss-repair fetch or accessed prefetch block. true if non-accessed prefetch block.
}**cache;

// Allocate memory for cache
struct cache_struct** allocate_cache(uint16_t CACHE_DEPTH, uint16_t N_WAY)
{
  struct cache_struct **created_array;

  // Allocate pointers in number equal to Cache depth to point to every row in the cache
  created_array = (struct cache_struct**) malloc(CACHE_DEPTH * sizeof(struct cache_struct*));

  uint64_t i;
  for (i = 0; i < CACHE_DEPTH; i++)
  {
    // Every row is allocated memory for cache lines equal to number of ways
    created_array[i] = (struct cache_struct*) malloc(N_WAY * sizeof(struct cache_struct));
  }

  return created_array;
}

// Initializa cache
void initialize_cache(void)
{
  uint16_t index, way;
  for(index=0; index<CACHE_DEPTH; index++)
  {
    for(way=0; way<N_WAY; way++)
    {
      cache[index][way].valid  = false;
      cache[index][way].dirty  = false;
      cache[index][way].prefetch_not_acc = false;
    }
  }
}

// Check if compare_tag matches with tag in any of the ways at the index. Return the way at which the compare_tag match is found.
bool check_hit_cache(uint64_t compare_tag, uint16_t index, uint16_t *way)
{
  uint16_t test_way;
  for(test_way=0; test_way<N_WAY; test_way++)
  {
    if((true == cache[index][test_way].valid) &&  (compare_tag == cache[index][test_way].tag))
    {
      *way = test_way;
      return true;
    }
  }

  return false;
}

// Returns true if cache is full else  provides the empty way
bool check_cache_full(uint16_t index, uint16_t *empty_way)
{
  uint16_t way;

  // Check if any way is empty
  for(way=0; way<N_WAY; way++)
  {
    if(false == cache[index][way].valid)
    {
      *empty_way = way;
      return false;
    }
  }

  return true; //Cache is full
}

// Returns true if cache is empty else  provides the first filled way
bool check_cache_empty(uint16_t index, uint16_t *filled_way)
{
  uint16_t way;

  // Check if any way is flled
  for(way=0; way<N_WAY; way++)
  {
    if(true == cache[index][way].valid)
    {
      *filled_way = way;
      return false;
    }
  }

  return true; //Cache is empty
}

// Find LRU way
uint16_t LRU_way_cache(uint16_t index)
{
  uint16_t way, lru_way;

  // Check LRU way
  for(way=0; way<N_WAY; way++)
  {
    if(0 == way)
    {
      lru_way = 0;
    }
    else if((true == cache[index][way].valid) && (cache[index][way].last_acc_time < cache[index][lru_way].last_acc_time))
    {
      lru_way = way;
    }
  }

  return lru_way;
}

// This structure templete represents 1-VC-line
struct VC_struct
{
  uint16_t cache_index; // Represents the index of the L1 cache to which this VC-line belongs
  uint64_t tag;
  bool valid;
  bool dirty;
  long inserted_time; // Indicates the time at which the VC-line was added to VC. This is needed to check for FI.
  bool prefetch_not_acc; // 'false' if miss-repair fetch or accessed prefetch block. true if non-accessed prefetch block.

  // VC-line represents the actual LRU cache-line for a given index.
}*VC;

// Allocate memory for VC
struct VC_struct* allocate_VC(uint8_t VC_DEPTH)
{
  struct VC_struct *created_array;

  // Allocate pointers in number equal to VC depth to point to every row in the VC
  created_array = (struct VC_struct*) malloc(VC_DEPTH * sizeof(struct VC_struct));

  return created_array;
}

// Initialize VC
void initialize_VC(void)
{
  uint8_t row;
  for(row=0; row<VC_DEPTH; row++)
  {
    VC[row].valid  = false;
    VC[row].dirty  = false;
    VC[row].prefetch_not_acc = false;
  }
}

// Check if compare_tag and cache_index matches with tag in any of the rows in VC. Return the row at which the compare_tag match is found.
bool check_hit_VC(uint64_t compare_tag, uint16_t cache_index, uint8_t *hit_row)
{
  uint8_t row;
  for(row=0; row<VC_DEPTH; row++)
  {
    //Check for cache_index and comapre_tag at VC[row] and VC[row] is valid
    if((true == VC[row].valid) && (VC[row].cache_index == cache_index) && (VC[row].tag == compare_tag))
    {
      *hit_row = row; // Need to return the row at which match was found
      return true; //cache_index and compare_tag were found
    }
  }

  return false; // either cache_index or compare_tag was not found in VC
}

// Returns true if VC is full else provides the empty row
bool check_VC_full(uint8_t *empty_row)
{
	uint8_t row;

  // Check if any row is empty
  for(row=0; row<VC_DEPTH; row++)
  {
    if(false == VC[row].valid)
    {
      *empty_row = row;
      return false;
    }
  }

	return true;
}

// Returns true if VC is empty else provides the filled row
bool check_VC_empty(uint8_t *filled_row)
{
	uint8_t row;

  // Check if any row is filled
  for(row=0; row<VC_DEPTH; row++)
  {
    if(true == VC[row].valid)
    {
      *filled_row = row;
      return false;
    }
  }

	return true; //VC is empty
}

// Get FI row in VC.
uint8_t first_in_row_VC(void)
{
  uint8_t row, first_in_row;

  for(row=0; row<VC_DEPTH; row++)
  {
    if(0 == row)
    {
      first_in_row = 0;
    }
    else if((true == VC[row].valid) && (VC[row].inserted_time < VC[first_in_row].inserted_time))
    {
      first_in_row = row;
    }
  }

  return first_in_row;
}


void strided_prefether(uint64_t block_address, struct cache_stats_t* p_stats)
{
  __int128_t current_stride = (__int128_t)block_address - (__int128_t)Last_Miss_Block_Addr;
  Last_Miss_Block_Addr = block_address;

  if(current_stride == pending_stride)
  {
    uint8_t i;
    for(i=1; i<=PREFETCH_DEGREE; i++)
    {
      __int128_t fetch_add = (__int128_t)block_address + i*(__int128_t)pending_stride;

      if( (fetch_add >= 0) && (fetch_add < (pow(2,64))) ) // fetch address is within 0 to (2^64 - 1)
      {
        prefetch((uint64_t)fetch_add, p_stats);

        // Prefetched 1-block. Update statistics
        p_stats->prefetched_blocks++;
        p_stats->bytes_transferred += pow(2,B);
      }
      else
      {
        break; // fetch address is not bound between  0 to (2^64 - 1)
      }
    }
  }

  pending_stride = current_stride;
}

void prefetch(uint64_t fetch_address, struct cache_stats_t* p_stats)
{
  uint64_t tag;
  uint16_t index, offset;
  Get_tag_index_offset(fetch_address, &tag, &index, &offset);

  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  long current_time;
  current_time = 1000*1000*1000 * (long)ts.tv_sec + (long)ts.tv_nsec;

  uint16_t way;
  uint8_t row;

  // Check if the address to be prefetched is already present in cache.
  if(true == check_hit_cache(tag, index, &way))
  {
    // Already present in Cache. Don't prefetch.
  }
  else if(true == check_hit_VC(tag, index, &row)) // Check if the address to be prefetched is already present in VC.
  {
    // Address to be prefetched is already present in VC. Move it to L1 as LRU
    // Copy VC row to temp
    uint64_t temp_tag = VC[row].tag;
    uint64_t temp_dirty = VC[row].dirty;

    // Address is stored in VC only if cache is full. Thus, no need to check if cache is not full.
    way = LRU_way_cache(index);

    // Put LRU way in cache to VC. Don't update FIFO stack
    VC[row].valid = true;
    VC[row].cache_index = index;
    VC[row].tag = cache[index][way].tag;
    VC[row].dirty = cache[index][way].dirty;
    //VC[row].inserted_time = current_time;
    VC[row].prefetch_not_acc = cache[index][way].prefetch_not_acc;

    // Copy VC_row from VC to Cache's lru way as LRU
    cache[index][way].valid = true;
    cache[index][way].tag = temp_tag;
    cache[index][way].dirty = temp_dirty;
    cache[index][way].prefetch_not_acc = true;
    // Do not update last_acc_time as prefetched block is LRU now

  }
  else if(false == check_cache_full(index, &way)) // Address not present in cache or VC. Check if cache is not full
  {
    // Cache is not full. Put prefetched address in empty way of cache
    cache[index][way].valid = true;
    cache[index][way].dirty = false;
    cache[index][way].tag = tag;
    cache[index][way].prefetch_not_acc = true;

    // Cache-way could be empty for prefetching as we can prefetch in multiple sets.
//    if(true == check_cache_empty(index, &way))
//    {
//      // Cache is empty. Update last accessed time as current time for way-0
//      cache[index][way].last_acc_time = current_time;
//    }
//    else
//    {
      // Cache is not empty. Prefetched block time is 1 less than LRU_way's time.
      uint16_t lru_way = LRU_way_cache(index);
      cache[index][way].last_acc_time = cache[index][lru_way].last_acc_time - 1;
//    }
  }
  else if(0 == VC_DEPTH) // Cache is full. Check if VC_DEPTH is 0
  {
    // VC is not there. Evict LRU from cache to make place for prefetched block
    way = LRU_way_cache(index);

    // If LRU is dirty update write-backs
    if(true == cache[index][way].dirty)
    {
      p_stats->write_backs++;

      // Update bytes transferred to memory. In writeback 2^B bytes are transferred to memory
      p_stats->bytes_transferred += pow(2,B);
    }

    // Put prefetched block in cache
    cache[index][way].valid = true;
    cache[index][way].dirty = false;
    cache[index][way].tag = tag;
    cache[index][way].prefetch_not_acc = true;
    // Do not update last_acc_time as prefetched block is LRU now
  }
  else // Cache is full. VC exists.
  {
    // Put prefetched address in LRU-way of cache. Put LRU-way of cache in empty or FI way of VC as last to be evicted

    if(true == check_VC_full(&row)) // Get empty row of VC if VC is not full.
    {
      // VC is full. Get FI row from VC
      row = first_in_row_VC();

      // If FI row is dirty, write-back to main memory
      if((true == VC[row].valid) && (true == VC[row].dirty))
      {
        p_stats->write_backs++;

        // Update bytes transferred to memory. In writeback 2^B bytes are transferred to memory
        p_stats->bytes_transferred += pow(2,B);
      }
    }

    // Get the LRU way from cache that will get moved to VC.
    way = LRU_way_cache(index);

    // LRU way from Cache stored in VC as FI. Update FIFO in VC
    VC[row].valid = true;
    VC[row].cache_index = index;
    VC[row].tag = cache[index][way].tag;
    VC[row].dirty = cache[index][way].dirty;
    VC[row].inserted_time = current_time;
    VC[row].prefetch_not_acc = cache[index][way].prefetch_not_acc;

    // Place prefetched tag to Cache's lru as lru
    cache[index][way].valid = true;
    cache[index][way].dirty = false;
    cache[index][way].tag = tag;
    cache[index][way].prefetch_not_acc = true;
    // Do not update last_acc_time as prefetched block is LRU now

  }

}

// Splits provided address into tag, index and offset
void Get_tag_index_offset(uint64_t address, uint64_t *tag, uint16_t *index, uint16_t *offset)
{
  uint64_t mask = 0;
  uint64_t i;
  uint64_t shift = 1;

  *tag = address >> (C-S);

  for(i=0; i<(C-S); i++)
  {
    mask ^= (shift<<i);
  }
  *index = ((address & mask) >> B);

  //*index = (address & uint64_t(pow(2,(C-S))-1) ) >> B;

  mask = 0;
  for(i=0; i<B; i++)
  {
    mask ^= (shift<<i);
  }
  *offset = address & mask;


  #ifdef debug
    printf("Address: %.16" PRIx64 ":: Tag: %.16" PRIx64 ":: Index: %.4" PRIx16 ":: Offset: %.4" PRIx16, address, *tag, *index, *offset);
  #endif

  #ifdef print_file
    fprintf(fw, "Address: %.16" PRIx64 ":: Tag: %.16" PRIx64 ":: Index: %.4" PRIx16 ":: Offset: %.4" PRIx16, address, *tag, *index, *offset);
  #endif

}


/**
 * Subroutine for initializing the cache. You many add and initialize any global or heap
 * variables as needed.
 * XXX: You're responsible for completing this routine
 *
 * @c The total number of bytes for data storage is 2^C
 * @b The size of a single cache line in bytes is 2^B
 * @s The number of blocks in each set is 2^S
 * @v The number of blocks in the victim cache is V
 * @k The prefetch distance is K
 */
void setup_cache(uint64_t c, uint64_t b, uint64_t s, uint64_t v, uint64_t k)
{
  C = c;
  B = b;
  S = s;

  CACHE_DEPTH = pow(2,(C-B-S));
  N_WAY = pow(2,S);

  VC_DEPTH = v;
  PREFETCH_DEGREE = k;

  #ifdef debug
    printf("Cache depth: %" PRIu16 "\nNo. of ways: %" PRIu16 "\nVC depth: %" PRIu16 "\nDegree of prefetch: %" PRIu16 "\n\n",
            CACHE_DEPTH, N_WAY, VC_DEPTH, PREFETCH_DEGREE);
  #endif

  #ifdef print_file
  	fprintf(fw, "Cache depth: %" PRIu16 "\nNo. of ways: %" PRIu16 "\nVC depth: %" PRIu16 "\nDegree of prefetch: %" PRIu16 "\n\n",
            CACHE_DEPTH, N_WAY, VC_DEPTH, PREFETCH_DEGREE);
  #endif

  cache = allocate_cache(CACHE_DEPTH, N_WAY);

  initialize_cache();

  if(0 != VC_DEPTH)
  {
  	VC = allocate_VC(VC_DEPTH);
  	initialize_VC();
	}

  Last_Miss_Block_Addr = 0;
  pending_stride = 0;

  #ifdef print_file
    fw = fopen("output_b.txt", "w");
    if (fw == NULL)
    {
      fprintf(stderr, "Can't open output file: Output.txt!\n");
      exit(1);
    }
  #endif

}


/**
 * Subroutine that simulates the cache one trace event at a time.
 * XXX: You're responsible for completing this routine
 *
 * @rw The type of event. Either READ or WRITE
 * @address  The target memory address
 * @p_stats Pointer to the statistics structure
 */
void cache_access(char rw, uint64_t address, struct cache_stats_t* p_stats)
{
  p_stats->accesses++;

  if('w' == rw)
  {
    p_stats->writes++;
  }
  else
  {
    p_stats->reads++;
  }

  uint64_t compare_tag;
  uint16_t index, offset, way;
  uint8_t row;
  // For given address get compare_tag, index and offset
  Get_tag_index_offset(address, &compare_tag, &index, &offset);

  struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
  long current_time;
  current_time = 1000*1000*1000 * (long)ts.tv_sec + (long)ts.tv_nsec;

  uint64_t block_address;
  block_address = address >> B;
  block_address <<= B;

  #ifdef debug
    printf(":: %lu", current_time);
  #endif

  #ifdef print_file
  	fprintf(fw, ":: Time: %lu", current_time);
  #endif


  // Check if compare tag is found at given index
  if(true == check_hit_cache(compare_tag, index, &way))
  {
    // Cache line was hit.

    // Check if this is a first time hit on prefetched block
    if(true == cache[index][way].prefetch_not_acc)
    {
      p_stats->useful_prefetches++;
      cache[index][way].prefetch_not_acc = false;
    }

    //update last accessed time to current time
    cache[index][way].last_acc_time = current_time;

    if('w' == rw)
    {

      #ifdef debug
        printf(":: WRITE HIT\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: WRITE HIT \n");
      #endif

      // Writing to the memory location. Since, we use WB scheme the cache line is dirty
      cache[index][way].dirty = true;
    }
    else
    {
      #ifdef debug
        printf(":: READ HIT\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: READ HIT \n");
      #endif
    }
  }
  else if(false == check_cache_full(index, &way)) //Miss on cache. Need to check if the cache is not full
  {
    // Cache is not full. Insert missed tag in empty way of cache
    p_stats->misses++;

    // For CPU only VC miss is an actual miss. Since, this is a miss for CPU increment VC misses.
    p_stats->vc_misses++;

    // We can insert the missed tag in the empty way
    cache[index][way].valid = true;
    cache[index][way].tag = compare_tag;
    cache[index][way].last_acc_time = current_time;
    cache[index][way].prefetch_not_acc = false;

    // Update bytes transferred from memory. In miss-repair2^B bytes are transferred from memory
    p_stats->bytes_transferred += pow(2,B);

    if('w' == rw) // Data is being written to cache
    {
      #ifdef debug
        //printf(":: WRITE MISS:: Put in Empty way of Cache\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: WRITE MISS:: Put in Empty way of Cache\n");
      #endif

      // New data being bought in for write hence it will be dirty
      cache[index][way].dirty = true;

      p_stats->write_misses++;
      p_stats->write_misses_combined++;
    }
    else // Data was being read from Cache
    {
      //New data being bought in for read hence not dirty
      cache[index][way].dirty = false;

      p_stats->read_misses++;
      p_stats->read_misses_combined++;
    }

    // Missed in cache. Call prefetcher
    strided_prefether(block_address, p_stats);

  }
  else if((0 == VC_DEPTH)) // No victim cache. Cache was missed and cache is full. Place tag in LRU way of cache
  {
    p_stats->misses++;

    // For CPU only VC miss is an actual miss. Since, this is a miss for CPU increment VC misses.
    p_stats->vc_misses++;

    // Check the LRU way to be replaced
    way = LRU_way_cache(index);

    // Block being flushed out of cache is dirty. Hence update write-backs
    if(true == cache[index][way].dirty)
    {
      p_stats->write_backs++;

      // Update bytes transferred to memory. In writeback 2^B bytes are transferred to memory
      p_stats->bytes_transferred += pow(2,B);
    }

    // New data is being bought to cache. Thus, it is valid with the new tag and current time
    cache[index][way].valid = true;
    cache[index][way].tag = compare_tag;
    cache[index][way].last_acc_time = current_time;
    cache[index][way].prefetch_not_acc = false;

    // Update bytes transferred from memory. In miss-repair2^B bytes are transferred from memory
    p_stats->bytes_transferred += pow(2,B);

    if('w' == rw)
    {
      #ifdef debug
        //printf(":: WRITE MISS:: Put in Empty LRU-way of Cache\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: WRITE MISS:: Put in Empty LRU-way of Cache\n");
      #endif

      //New data being bought in for write hence it will be dirty
      cache[index][way].dirty = true;

      p_stats->write_misses++;
      p_stats->write_misses_combined++;
    }
    else //New data being bought in for read hence not dirty
    {
      #ifdef debug
        //printf(":: READ MISS:: Put in Empty LRU-way of Cache\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: READ MISS:: Put in Empty LRU-way of Cache\n");
      #endif

      cache[index][way].dirty = false;

      p_stats->read_misses++;
      p_stats->read_misses_combined++;
    }

    // Missed in cache. Call prefetcher
    strided_prefether(block_address, p_stats);

  }
  else if(true == check_hit_VC(compare_tag, index, &row)) //Cache is full. VC exists. Check in VC .
  {
    // Hit in VC. Need to to swap between VC and cache LRU way
    p_stats->misses++;

    // Check if this is a first time hit on prefetched block
    if(true == VC[row].prefetch_not_acc)
    {
      p_stats->useful_prefetches++;
      VC[row].prefetch_not_acc = false;
    }

    // Copy VC[row] into temporary storage
    uint64_t temp_tag = VC[row].tag;
    uint64_t temp_dirty = VC[row].dirty;

    // Get the LRU way from cache that will get swaped.
    way = LRU_way_cache(index);

    // LRU way from Cache stored in VC as FI
    VC[row].valid = true;
    VC[row].cache_index = index;
    VC[row].tag = cache[index][way].tag;
    VC[row].dirty = cache[index][way].dirty;
    //VC[row].inserted_time = current_time;
    VC[row].prefetch_not_acc = cache[index][way].prefetch_not_acc;

    // Copy VC_row from VC to Cache's lru way as MRU
    cache[index][way].valid = true;
    cache[index][way].tag = temp_tag;
    cache[index][way].dirty = temp_dirty;
    cache[index][way].last_acc_time = current_time;
    cache[index][way].prefetch_not_acc = false;

    if('w' == rw) // Data is being written
    {
      #ifdef debug
        printf(" :: WRITE Hit in VC\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: WRITE Hit in VC\n");
      #endif

      // New data being bought in for write hence it will be dirty
      cache[index][way].dirty = true;

      p_stats->write_misses++;
    }
    else // Data was being read from Cache
    {
      #ifdef debug
        printf(" :: READ Hit in VC\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: READ Hit in VC\n");
      #endif

      p_stats->read_misses++;
    }

    // Missed in cache. Call prefetcher
    strided_prefether(block_address, p_stats);

  }
  else if(false == check_VC_full(&row)) // Miss in VC. Check if VC is not full.
  {
    // VC is not full. Need to get LRU way from cache to VC and missed tag in LRU way of cache
    p_stats->misses++;
    p_stats->vc_misses++;

    // Get the LRU way from cache that will get moved to VC.
    way = LRU_way_cache(index);

    // LRU way from Cache stored in VC as FI
    VC[row].valid = true;
    VC[row].cache_index = index;
    VC[row].tag = cache[index][way].tag;
    VC[row].dirty = cache[index][way].dirty;
    VC[row].inserted_time = current_time;
    VC[row].prefetch_not_acc = cache[index][way].prefetch_not_acc;

    // Place missed tag to Cache's lru way as MRU
    cache[index][way].valid = true;
    cache[index][way].tag = compare_tag;
    cache[index][way].last_acc_time = current_time;
    cache[index][way].prefetch_not_acc = false;

    // Update bytes transferred from memory. In miss-repair2^B bytes are transferred from memory
    p_stats->bytes_transferred += pow(2,B);

    if('w' == rw) // Data is being written
    {
      #ifdef debug
        //printf(" :: WRITE MISS in VC\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: WRITE MISS in VC\n");
      #endif

      // New data being bought in for write hence it will be dirty
      cache[index][way].dirty = true;

      p_stats->write_misses++;
      p_stats->write_misses_combined++;
    }
    else // Data was being read from Cache
    {
      #ifdef debug
        //printf(" :: READ MISS in VC\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: READ MISS in VC\n");
      #endif

      // New data being bought in for read hence it will not be dirty
      cache[index][way].dirty = false;

      p_stats->read_misses++;
      p_stats->read_misses_combined++;
    }

    // Missed in cache. Call prefetcher
    strided_prefether(block_address, p_stats);

  }
  else //Miss in VC and VC is full.
  {
    // Flush FI row in VC out. Move LRU-way from cache to VC. Get missed tag in LRU-way of cache
    p_stats->misses++;
    p_stats->vc_misses++;

    // Get FI row from VC
    row = first_in_row_VC();

    // If FI row is dirty, write-back to main memory
    if(true == VC[row].dirty)
    {
      p_stats->write_backs++;

      // Update bytes transferred to memory. In writeback 2^B bytes are transferred to memory
      p_stats->bytes_transferred += pow(2,B);
    }

    // Get the LRU way from cache that will get moved to VC.
    way = LRU_way_cache(index);

    // Move LRU-way from cache to FI row VC. Update FIFO
    VC[row].valid = true;
    VC[row].cache_index = index;
    VC[row].tag = cache[index][way].tag;
    VC[row].dirty = cache[index][way].dirty;
    VC[row].inserted_time = current_time;
    VC[row].prefetch_not_acc = cache[index][way].prefetch_not_acc;

    // Place missed tag to Cache's lru way as MRU
    cache[index][way].valid = true;
    cache[index][way].tag = compare_tag;
    cache[index][way].last_acc_time = current_time;
    cache[index][way].prefetch_not_acc = false;

    // Update bytes transferred from memory. In miss-repair2^B bytes are transferred from memory
    p_stats->bytes_transferred += pow(2,B);

    if('w' == rw) // Data is being written
    {
      #ifdef debug
        //printf(" :: WRITE MISS in VC\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: WRITE MISS in VC\n");
      #endif

      // New data being bought in for write hence it will be dirty
      cache[index][way].dirty = true;

      p_stats->write_misses++;
      p_stats->write_misses_combined++;
    }
    else // Data was being read from Cache
    {
      #ifdef debug
        //printf(" :: READ MISS in VC\n");
      #endif

      #ifdef print_file
        fprintf(fw, ":: READ MISS in VC\n");
      #endif

      // New data being bought in for read hence it will not be dirty
      cache[index][way].dirty = false;

      p_stats->read_misses++;
      p_stats->read_misses_combined++;
    }

    // Missed in cache. Call prefetcher
    strided_prefether(block_address, p_stats);

  }

}



/**
 * Subroutine for cleaning up any outstanding memory operations and calculating overall statistics
 * such as miss rate or average access time.
 * XXX: You're responsible for completing this routine
 *
 * @p_stats Pointer to the statistics structure
 */
void complete_cache(struct cache_stats_t *p_stats)
{
  #ifdef debug
    printf("\n");
  #endif

  #ifdef print_file
    fclose(fw);
  #endif

  p_stats->hit_time        = 2 + (double)0.2*S;
  p_stats->miss_penalty    = 200;
  p_stats->miss_rate       = ((double)p_stats->misses/(double)p_stats->accesses);
  p_stats->avg_access_time = p_stats->hit_time + ((double)p_stats->vc_misses/(double)p_stats->accesses) * p_stats->miss_penalty;

  uint16_t i;
  //Freeing up allocated memory
  for(i=0; i<CACHE_DEPTH; i++)
  {
    free(cache[i]);
  }

  free(cache);
  free(VC);

}
