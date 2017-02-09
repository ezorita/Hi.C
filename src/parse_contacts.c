#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <search.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <libgen.h>
#include <ctype.h>

#define    HIC_FORMAT 1
#define COOLER_FORMAT 2

#define FORMAT HIC_FORMAT

#define MIN_SCORE 10
#define HASH_SIZE 2048
#define MAX_INSERT_SIZE 2000
#define MIN_MAPQ 20
#define MAX_OVERLAP 4

#define FLAG_MULTISEGMENT   0x001
#define FLAG_PROPALIGN      0x002
#define FLAG_UNMAPPED       0x004
#define FLAG_NEXT_UNMAPPED  0x008
#define FLAG_REVCOMP        0x010
#define FLAG_NEXT_REVCOMP   0x020
#define FLAG_FORWARD_READ   0x040
#define FLAG_REVERSE_READ   0x080
#define FLAG_SECONDARY      0x100
#define FLAG_FILTERED       0x200
#define FLAG_PCR_DUPLICATE  0x400
#define FLAG_SUPPL_ALIGN    0x800

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

#define lowercase(s) for(char * p = s;*p;++p) *p=tolower(*p)

// Struct definitions.

typedef struct {
   char * seqname;
   char * chr;
   char * cigar;
   short  flag;
   short  mapq;
   int    score;
   long   locus;
} sam_t;

typedef struct {
   char * chr;
   long   beg_ref;
   long   end_ref;
   long   beg_frag;
   long   end_frag;
   int    beg_read;
   int    end_read;
   int    frag_id;
   int    rc;
   int    mapq;
} map_t;

typedef struct {
   int      pos;
   int      max; 
   sam_t ** buf;
} samstack_t;

typedef struct {
   int   pos;
   int   max;
   map_t map[];
} mapstack_t;

typedef struct {
   int beg_clip;
   int end_clip;
   int matches;
   int insertions;
   int deletions;
} cigar_t;

typedef struct {
   char * chr;
   int    cnt;
   int  * re_site;
} isdchr_t;



// Function headers.
void           parse_sam    (sam_t * sam, char * samline);
cigar_t        parse_cigar  (char *);
int            parse_contact(samstack_t *, struct hsearch_data *, int, int);
sam_t        * new_sam      (void);
samstack_t   * new_samstack (int);
mapstack_t   * new_mapstack (int);
int            sam_push     (sam_t *, samstack_t *);
int             map_push     (map_t, mapstack_t **);
void           sam_destroy  (sam_t * sam);
int            find_pe_contacts (mapstack_t  * fw, mapstack_t  * rv, mapstack_t ** dst);
void           place_in_read (mapstack_t * src, mapstack_t * dst);

// Restriction enzyme functions.
int      parse_isd    (int re_fd, char* re_name, struct hsearch_data* htable);
int      bisection    (int* data, int beg, int end, int target);
int      read_enzyme_db (char *, char *, struct hsearch_data *);
void     fill_re_fragment_info (map_t *, struct hsearch_data *);

// Sort compar functions.
int            sam_by_score_desc (const void *, const void *);
int            map_by_read_beg   (const void *, const void *);

// Stats variables.
long valid         = 0;
long single_read   = 0;
long unmapped      = 0;
long repeats       = 0;
long self_filter   = 0;
long dangling      = 0;
long unknown       = 0;
long insert_filter = 0;


int main(int argc, char *argv[])
{
   // Parse params.
   if (argc < 4) {
      fprintf(stderr, "usage: %s <organism> <RE> <hic-pe.sam> [mapq >= 20] [ins_size <= 2000]\n", argv[0]);
      exit(1);
   }

   fprintf(stderr, "open files...");

   // Get args.
   char * organism = argv[1];
   char * re_name  = argv[2];
   int    min_mapq = MIN_MAPQ;
   int    max_insz = MAX_INSERT_SIZE;
   if (argc > 4) min_mapq = atoi(argv[4]);
   if (argc > 5) max_insz = atoi(argv[5]);

   // Open files.
   FILE * fin = fopen(argv[3], "r");
   if (fin == NULL) {
      fprintf(stderr, "error opening file: %s\n", argv[3]);
      exit(1);
   }

   // Read database.
   fprintf(stderr, "ok\nloading RE database...");
   struct hsearch_data htable = {0};
   hcreate_r(HASH_SIZE, &htable);
   read_enzyme_db(organism, re_name, &htable);
   fprintf(stderr, "ok\nparsing sam file...");

   // Read lines.
   size_t lines = 0;
   size_t bufsize = 200;
   char * line = malloc(bufsize);
   ssize_t bytes = 0;
   char sam_last;
   // SAM buffers.
   samstack_t * stack = new_samstack(100);
   sam_t      * sam   = new_sam();
   

   // Skip header until first read.
   do {
      bytes = getline(&line, &bufsize, fin);
   } while (bytes > 0 && line[0] == '@');
   // End of file.
   if (bytes < 0) {
      fprintf(stderr,"error: input file is empty.\n");
      exit(1);
   }
   // Parse line.
   parse_sam(sam, line);

   // File loop.
   char * seqname = malloc(1000);
   do {
      // 0. Reset variables.
      stack->pos = 0;
      // 1. Save reference read.
      strcpy(seqname, sam->seqname);
      sam_push(sam, stack);
      // 2. Store all reads with the same seqname.
      do {
         bytes = getline(&line, &bufsize, fin);
         // End of file.
         if (bytes < 0) break;
         // Parse line.
         parse_sam(sam, line);
         // Compare seqnames.
         if (strcmp(seqname, sam->seqname) != 0)
            break;
         sam_push(sam, stack);
      } while (1);
      // 3. Process reads.
      parse_contact(stack, &htable, min_mapq, max_insz);
   } while (bytes > 0);

   fprintf(stderr, "ok\n\nValid pairs:            \t%ld\n", valid);
   fprintf(stderr, "Invalid pairs:          \t%ld\n", 
           insert_filter+self_filter+single_read+unmapped+repeats+dangling+unknown);
   fprintf(stderr, " - Repeats:             \t%ld\n", repeats);
   fprintf(stderr, " - Dangling ends:       \t%ld\n", dangling);
   fprintf(stderr, " - Self ligated:        \t%ld\n", self_filter);
   fprintf(stderr, " - One read mapped:     \t%ld\n", single_read);
   fprintf(stderr, " - Unmapped:            \t%ld\n", unmapped);
   fprintf(stderr, " - Insert size (>%dbp):\t%ld\n", max_insz, insert_filter);
   fprintf(stderr, " - Unknown event:       \t%ld\n", unknown);
   free(stack);
   free(sam);

   return 0;
}
 
int
parse_contact
(
 samstack_t * stack,
 struct hsearch_data * htable,
 int min_mapq,
 int max_insz
)
{
   // Sort sam by score.
   qsort(stack->buf, stack->pos, sizeof(sam_t*), sam_by_score_desc);

   // Note: Contacts in the same reads are directly accepted (if Q > thr).
   // Contacts between reads must satisfy insert size restrictions as well.
   mapstack_t * mapf = new_mapstack(10);
   mapstack_t * mapr = new_mapstack(10);
   for (int i = 0; i < stack->pos; i++) {
      sam_t * sam = stack->buf[i];
      int f = sam->flag;
      // Filter unmapped reads.
      if ((f & FLAG_UNMAPPED) || !(f & FLAG_MULTISEGMENT))
         continue;

      // Define map coordinates 5'->3'.
      cigar_t cigar = parse_cigar(sam->cigar);
      int   revcomp = sam->flag & FLAG_REVCOMP;
      map_t map = (map_t) {
         .chr      = sam->chr,
         .beg_ref  = sam->locus,
         .end_ref  = sam->locus + cigar.matches + cigar.deletions,
         .beg_frag = -1,
         .end_frag = -1,
         .beg_read = (revcomp ? cigar.end_clip : cigar.beg_clip),
         .end_read = (revcomp ? cigar.end_clip : cigar.beg_clip) + cigar.matches + cigar.insertions,
         .frag_id  = -1,
         .rc       = revcomp, 
         .mapq     = sam->mapq
      };

      // Classify forward and reverse reads.
      if (f & FLAG_FORWARD_READ)
         map_push(map, &mapf);
      else 
         map_push(map, &mapr);
   }

   mapstack_t * mf = new_mapstack(max(1,mapf->pos));
   mapstack_t * mr = new_mapstack(max(1,mapr->pos));

   // Apply filters.
   if (mapr->pos + mapf->pos == 0) {
      unmapped++;
      goto free_and_return;
   }
   else if (mapr->pos + mapf->pos == 1) {
      single_read++;
      goto free_and_return;
   }


   // Place mappings in read. (First to be placed is best alignment score, then the others if they fit).
   place_in_read(mapf, mf);
   place_in_read(mapr, mr);


   // Filter by quality and fill restriction enzyme fragment info.
   mapf->pos = mapr->pos = 0;
   for (int i = 0; i < mf->pos; i++) {
      if (mf->map[i].mapq >= min_mapq) {
         mapf->map[mapf->pos] = mf->map[i];
         fill_re_fragment_info(mapf->map+(mapf->pos++), htable);
      }
   }

   for (int i = 0; i < mr->pos; i++) {
      if (mr->map[i].mapq >= min_mapq) {
         mapr->map[mapr->pos] = mr->map[i];
         fill_re_fragment_info(mapr->map+(mapr->pos++), htable);
      }
   }

   if (mapr->pos + mapf->pos < 2) {
      repeats++;
      goto free_and_return;
   }

   // Millor primer: Join fragments between reads. (aixo eliminara duplicats)
   // Despres: Imprimir els contactes tal qual.
   int insert_size = find_pe_contacts(mapf, mapr, &mf);
   if (insert_size > max_insz) {
      insert_filter++;
      goto free_and_return;
   }
   
   // Output contacts.
   for (int i = 0; i < mf->pos-1; i++) {
      for (int j = i+1; j < mf->pos; j++) {
         valid++;
         map_t m1 = mf->map[i];
         map_t m2 = mf->map[j];
         int chrcmp = strcmp(m1.chr, m2.chr);
         if (chrcmp > 0 || (chrcmp == 0 && m2.beg_ref < m1.beg_ref)) {
            m2 = mf->map[i];
            m1 = mf->map[j];
         }
#if FORMAT == HIC_FORMAT
         fprintf(stdout, "%s %d %s %ld %d %d %s %ld %d %d %d\n",
                 stack->buf[0]->seqname,
                 m1.rc ? 1 : 0,
                 m1.chr,
                 m1.beg_ref,
                 m1.frag_id,
                 m2.rc ? 1 : 0,
                 m2.chr,
                 m2.beg_ref,
                 m2.frag_id,
                 m1.mapq,
                 m2.mapq
                 );
#elif FORMAT == COOLER_FORMAT
         fprintf(stdout, "%s\t%ld\t%s\t%s\t%ld\t%s\n",
                 m1.chr,
                 m1.beg_ref,
                 m1.rc ? "-" : "+",
                 m2.chr,
                 m2.beg_ref,
                 m2.rc ? "-" : "+"
                 );
#else
         fprintf(stdout, "%s\t%s\t%ld\t%d\t%ld\t%ld\t%ld\t%s\t%ld\t%d\t%ld\t%ld\t%ld\n",
                 stack->buf[0]->seqname,
                 m1.chr,
                 m1.beg_ref,
                 m1.rc ? 1 : 0,
                 m1.end_ref - m1.beg_ref,
                 m1.beg_frag,
                 m1.end_frag,
                 m2.chr,
                 m2.beg_ref,
                 m2.rc ? 1 : 0,
                 m2.end_ref - m2.beg_ref,
                 m2.beg_frag,
                 m2.end_frag
                 );
#endif

      }
   }

 free_and_return:
   free(mapf);
   free(mapr);
   free(mr);
   free(mf);
   return 0;
}

int
find_pe_contacts
(
 mapstack_t  * fw,
 mapstack_t  * rv,
 mapstack_t ** dst
)
{
   int insert_size = 0;
   (*dst)->pos = 0;
   // Sort by position in read.
   qsort(fw->map, fw->pos, sizeof(map_t), map_by_read_beg);
   qsort(rv->map, rv->pos, sizeof(map_t), map_by_read_beg);

   int inner_merged = 0;
   int outer_merged = 0;
   int self_ligation = 0;
   int unknown_event = 0;
   if (fw->pos && rv->pos) {
      int multi_map = (fw->pos > 1) && (rv->pos > 1);
      int check_outer_loop = 0;
      // Check inner loop.
      map_t ifw = fw->map[fw->pos-1];
      map_t irv = rv->map[rv->pos-1];
      if (ifw.frag_id == irv.frag_id) {
         if (ifw.rc != irv.rc) {
            // Check whether fragment is contiguous or self-ligated.
            if ((ifw.rc && (ifw.end_ref < irv.beg_ref)) || (irv.rc && (ifw.beg_ref > irv.end_ref)))
               self_ligation = 1;
            // Merge fragments.
            ifw.beg_ref = min(ifw.beg_ref, irv.beg_ref);
            ifw.end_ref = max(ifw.end_ref, irv.end_ref);
            ifw.mapq    = max(ifw.mapq, irv.mapq);
            insert_size = ifw.end_ref - ifw.beg_ref + 1;
            map_push(ifw, dst);
            inner_merged = 1;
         } else {
            // What is this? Same fragment sequenced in the same direction?
            // This would require two exactly equal molecules or a broken
            // molecule that flipped and ligated to itself --> classify as 
            // 'unknown'.
            unknown_event = 1;
         }
      } else {
         check_outer_loop = 1;
         // Compute insert size by sum of fragments.
         if (ifw.rc)
            insert_size += ifw.end_frag - ifw.beg_ref;
         else
            insert_size += ifw.end_ref - ifw.beg_frag;
         if (irv.rc)
            insert_size += irv.end_ref - irv.beg_frag;
         else
            insert_size += irv.end_frag - irv.beg_ref;
      }
      if (check_outer_loop || multi_map) {
         // Check molecule outer loop.
         map_t ofw = fw->map[0];
         map_t orv = rv->map[0];
         if (ofw.frag_id == orv.frag_id) {
            if (ofw.rc != orv.rc) {
               ofw.beg_ref = min(ofw.beg_ref, orv.beg_ref);
               ofw.end_ref = max(ofw.end_ref, orv.end_ref);
               ofw.mapq    = max(ofw.mapq, orv.mapq);
               map_push(ofw, dst);
            } else {
               // Outer loop does not make sense. What do we do?
               // keep the one with greater mapq.
               if (ofw.mapq < orv.mapq)
                  map_push(orv, dst);
               else
                  map_push(ofw, dst);
            }
            outer_merged = 1;
         }
      }
   }
   // Add other fragments.
   mapstack_t * tmp = new_mapstack(fw->pos+rv->pos);
   int beg = (outer_merged ? 1 : 0);
   int end_fw = fw->pos - (inner_merged ? 1 : 0);
   int end_rv = rv->pos - (inner_merged ? 1 : 0);
   // Fill tmp with remaining fragments.
   for (int i = beg ; i < end_fw ; i++)
      tmp->map[tmp->pos++] = fw->map[i];
   for (int i = beg ; i < end_rv ; i++) {
      rv->map[i].rc = !rv->map[i].rc;
      tmp->map[tmp->pos++] = rv->map[i];
   }
   // Add fragments without repetition.
   for (int i = 0; i < tmp->pos; i++) {
      int add = 1;
      for (int j = i + 1; j < tmp->pos; j++) {
         if (tmp->map[i].frag_id == tmp->map[j].frag_id) {
            if (tmp->map[i].mapq > tmp->map[j].mapq)
               tmp->map[j] = tmp->map[i];
            add = 0;
            break;
         }
      }
      if (add)
         map_push(tmp->map[i], dst);
   }

   if ((*dst)->pos == 1) {
      if (self_ligation)
         self_filter++;
      else if (unknown_event)
         unknown++;
      else if (inner_merged)
         dangling++;
      else
         single_read++;
   }

   free(tmp);

   return insert_size;
}

void
place_in_read
(
 mapstack_t * src,
 mapstack_t * dst
)
{
   if (src->pos == 0) return;
   // Add first.
   dst->map[dst->pos++] = src->map[0];
   for (int i = 1; i < src->pos; i++) {
      map_t new = src->map[i];
      int insert = 1;
      for (int j = 0; j < dst->pos; j++) {
         map_t old = dst->map[j];
         // Compute overlap.
         int beg = max(old.beg_read, new.beg_read);
         int end = min(old.end_read, new.end_read);
         int overlap = max(0,end-beg+1);         
         if (overlap > MAX_OVERLAP) {
            insert = 0;
            break;
         }
      }
      if (insert)
         dst->map[dst->pos++] = new;
   }
}

void
fill_re_fragment_info
(
 map_t * map,
 struct hsearch_data * htable
)
{
   // Get chromosome RE sites.
   ENTRY * item;
   hsearch_r((ENTRY){.key = map->chr}, FIND, &item, htable);
   
   if (item == NULL) {
      fprintf(stderr, "warning: chromosome not found in digestion file: %s. RE info set to -1.\n", map->chr);
      return;
   }
   isdchr_t * ref = (isdchr_t *) item->data;
   // Find fragment by bisection.
   int idx = bisection(ref->re_site, 0, ref->cnt-1, map->beg_ref);
   map->beg_frag = ref->re_site[idx];
   map->end_frag = ref->re_site[idx+1];
   map->frag_id  = idx;
   
   return;
}

int
bisection
(
 int * data,
 int   beg,
 int   end,
 int   target
 )
{
   if (end - beg < 2) return beg;
   int mid = (beg+end)/2;
   if (target < data[mid]) end = mid;
   else if (target > data[mid]) beg = mid;
   else return mid;
   
   return bisection(data,beg,end,target);
}

cigar_t
parse_cigar
(
 char * str
)
{
   cigar_t cigar = (cigar_t){.beg_clip = 0, .end_clip = 0, .matches = 0, .insertions = 0, .deletions = 0};

   if (str[0] == '*') return cigar;

   char * num = str;
   int len = strlen(str);
   for (int i = 0; i < len; i++) {
      if (str[i] == 'H' || str[i] == 'S') {
         str[i] = 0;
         if (cigar.matches + cigar.insertions + cigar.deletions)
            cigar.end_clip += atoi(num);
         else
            cigar.beg_clip += atoi(num);
         num = str+i+1;
      } else if (str[i] == 'M') {
         str[i] = 0;
         cigar.matches += atoi(num);
         num = str+i+1;
      } else if (str[i] == 'I') {
         str[i] = 0;
         cigar.insertions += atoi(num);
         num = str+i+1;
      } else if (str[i] == 'D') {
         str[i] = 0;
         cigar.deletions += atoi(num);
         num = str+i+1;
      }
   }
   return cigar;
}

void
parse_sam
(
 sam_t * sam,
 char  * samline
)
{
   strcpy(sam->seqname,strtok(samline, "\t"));
   sam->flag = atoi(strtok(NULL,"\t"));
   strcpy(sam->chr, strtok(NULL,"\t"));
   sam->locus = atoi(strtok(NULL,"\t"));
   sam->mapq = atoi(strtok(NULL,"\t"));
   strcpy(sam->cigar, strtok(NULL,"\t"));
   char * str;
   while (str = strtok(NULL,"\t")) {
      if (str[0] == 'A' && str[1] == 'S') {
         sam->score = atoi(str+5);
         break;
      }
   }
}

void
samcpy
(
 sam_t * sam_d,
 sam_t * sam_s
)
{
   strcpy(sam_d->seqname,sam_s->seqname);
   sam_d->flag = sam_s->flag;
   strcpy(sam_d->chr, sam_s->chr);
   sam_d->locus = sam_s->locus;
   sam_d->mapq = sam_s->mapq;
   sam_d->score = sam_s->score;
   strcpy(sam_d->cigar, sam_s->cigar);
}


sam_t *
new_sam
(
 void
)
{
   sam_t * sam = malloc(sizeof(sam_t));
   if (!sam) return NULL;
   sam->seqname = malloc(100);
   sam->flag = 0;
   sam->chr = malloc(100);
   sam->locus = 0;
   sam->mapq = 0;
   sam->cigar = malloc(100);

   return sam;
}

void
sam_destroy
(
 sam_t * sam
)
{
   free(sam->seqname);
   free(sam->chr);
   free(sam->cigar);
   free(sam);
}

samstack_t *
new_samstack
(
 int size
)
{
   samstack_t * stack = malloc(sizeof(samstack_t));
   sam_t     ** buf   = malloc(size*sizeof(sam_t*));
   if (stack == NULL || buf == NULL)
      return NULL;

   for (int i = 0; i < size; i++) {
      buf[i] = new_sam();
   }

   stack->max = size;
   stack->pos = 0;
   stack->buf = buf;

   return stack;
}

void
samstack_destroy
(
 samstack_t * stack
)
{
   for (int i = 0; i < stack->max; i++) 
      sam_destroy(stack->buf[i]);

   free(stack);
}

mapstack_t *
new_mapstack
(
 int size
)
{
   mapstack_t * stack = malloc(sizeof(mapstack_t) + size * sizeof(map_t));
   if (stack == NULL)
      return NULL;
   stack->max = size;
   stack->pos = 0;

   return stack;
}

int
sam_push
(
 sam_t      * sam,
 samstack_t * stack
)
{
   if (stack->pos >= stack->max) {
      int newsize = 2*stack->max;
      stack->buf= realloc(stack->buf, newsize*sizeof(sam_t *));
      if (stack->buf == NULL)
         return 1;
      stack->max = newsize;
   }
   samcpy(stack->buf[stack->pos++], sam);
   return 0;
}

int
map_push
(
 map_t         map,
 mapstack_t ** stackp
)
{
   mapstack_t * stack = *stackp;

   if (stack->pos >= stack->max) {
      int newsize = 2*stack->max;
      stack = *stackp = realloc(stack, sizeof(mapstack_t) + newsize * sizeof(map_t));
      if (stack == NULL)
         return 1;
      stack->max = newsize;
   }
   stack->map[stack->pos++] = map;
   return 0;
}

int 
sam_by_score_desc
(
 const void * ap,
 const void * bp
)
{
   sam_t * a = *((sam_t **) ap);
   sam_t * b = *((sam_t **) bp);
   
   if (a->score <= b->score) return 1;
   else return -1;
}

int 
map_by_read_beg
(
 const void * ap,
 const void * bp
)
{
   map_t * a = (map_t *) ap;
   map_t * b = (map_t *) bp;
   
   if (a->beg_read >= b->beg_read) return 1;
   else return -1;
}

int 
map_by_read_end
(
 const void * ap,
 const void * bp
)
{
   map_t * a = (map_t *) ap;
   map_t * b = (map_t *) bp;
   
   if (a->end_read >= b->end_read) return 1;
   else return -1;
}

int
read_enzyme_db
(
 char                * organism,
 char                * re_name,
 struct hsearch_data * htable
)
{
   // Open digest files.
   // To lowercase.
   lowercase(re_name);
   char * db_path = malloc(strlen(re_name)+strlen(organism)+9);
   sprintf(db_path, "db/%s/%s.isd", organism, re_name);

   // Open file and store fd.
   int re_fd = open(db_path, O_RDONLY);
   free(db_path);
   if (re_fd < 0) {
      fprintf(stderr, "error while opening RE database: %s.\n",db_path);
      exit(1);
   }

   // Parse digest files.
   parse_isd(re_fd, re_name, htable);

   close(re_fd);
}  

int
parse_isd
(
 int         re_fd,
 char      * re_name,
 struct hsearch_data * htable
)
{
   struct stat sb;
   char * p, * pmap;
   
   // mmap file.
   if (fstat(re_fd, &sb) == -1) {
      fprintf(stderr, "error reading digestion file (fstat).\n");
      exit(1);
   }

   pmap = mmap(0, sb.st_size, PROT_READ, MAP_SHARED, re_fd, 0);
   if (pmap == NULL) {
      fprintf(stderr, "error reading digestion file (mmap).\n");
      exit(1);
   }

   p = pmap;
   // Read RE information.
   p += strlen(p)+1;
   p += sizeof(int);
   p += sizeof(int);

   // Get number of chromosomes.
   int nchrom = *((int *)p);
   p += sizeof(int);
   
   // Parse each chromosome.
   int re_name_len = strlen(re_name);
   for (int i = 0; i < nchrom; i++) {
      // Create isd chromosome entry.
      isdchr_t * isdchr = malloc(sizeof(isdchr_t));

      // Insert isdchr in hash table (key is chromosome name).
      ENTRY * item;
      hsearch_r((ENTRY){.key = p, .data = isdchr}, ENTER, &item, htable);
     
      // Chromosome name.
      isdchr->chr = p;
      p += strlen(isdchr->chr)+1;
      // Number of RE sites.
      isdchr->cnt = *((int *)p);
      p += sizeof(int);
      // RE site list.
      isdchr->re_site = ((int *)p);
      p += isdchr->cnt*sizeof(int);
   }

   return 0;
}
      
