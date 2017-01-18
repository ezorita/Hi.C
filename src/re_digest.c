#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <search.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <libgen.h>
#include <ctype.h>

#define HASH_SIZE 2048
#define MAX_FRAGMENT_SIZE 2000

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

#define lowercase(s) for(char * p = s;*p;++p) *p=tolower(*p)

#define HELP_MSG "Need help? Call 911 modafacka.\n"


/* Nucleic acid notation:
** 
** A [0b0001]: Adenine
** C [0b0010]: Cytosine
** G [0b0100]: Guanine
** T [0b1000]: Thymine
** W [0b1001]: Weak (A or T)
** S [0b0110]: Strong (C or G)
** M [0b0011]: Amino (A or C)
** K [0b1100]: Keto (G or T)
** R [0b0101]: Purine (A or G)
** Y [0b1010]: Pyridimine (C or T)
** B [0b1110]: not A
** D [0b1101]: not C
** H [0b1011]: not G
** V [0b0111]: not T
** N [0b1111]: Any nucleotide
*/

#define NT_A 0b0001
#define NT_C 0b0010
#define NT_G 0b0100
#define NT_T 0b1000
#define NT_W 0b1001
#define NT_S 0b0110
#define NT_M 0b0011
#define NT_K 0b1100
#define NT_R 0b0101
#define NT_Y 0b1010
#define NT_B 0b1110
#define NT_D 0b1101
#define NT_H 0b1011
#define NT_V 0b0111
#define NT_N 0b1111
#define NT_Z 0b0000

// Reverse complements

#define RC_A NT_T
#define RC_C NT_G
#define RC_G NT_C
#define RC_T NT_A
#define RC_W NT_W
#define RC_S NT_S
#define RC_M NT_K
#define RC_K NT_M
#define RC_R NT_Y
#define RC_Y NT_R
#define RC_B NT_V
#define RC_D NT_H
#define RC_H NT_D
#define RC_V NT_B
#define RC_N NT_N


// Variable definitions.

const char re_nt[128] = {
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,NT_A,NT_B,NT_C,NT_D,0,
   0,NT_G,NT_H,0,0,NT_K,0,NT_M,NT_N,0,
   0,0,NT_R,NT_S,NT_T,0,NT_V,NT_W,0,NT_Y,
   0,0,0,0,0,0,NT_A,NT_B,NT_C,NT_D,0,
   0,NT_G,NT_H,0,0,NT_K,0,NT_M,NT_N,0,
   0,0,NT_R,NT_S,NT_T,0,NT_V,NT_W,0,NT_Y,
   0,0,0,0,0
};

const char re_rc[128] = {
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,RC_A,RC_B,RC_C,RC_D,0,
   0,RC_G,RC_H,0,0,RC_K,0,RC_M,RC_N,0,
   0,0,RC_R,RC_S,RC_T,0,RC_V,RC_W,0,RC_Y,
   0,0,0,0,0,0,RC_A,RC_B,RC_C,RC_D,0,
   0,RC_G,RC_H,0,0,RC_K,0,RC_M,RC_N,0,
   0,0,RC_R,RC_S,RC_T,0,RC_V,RC_W,0,RC_Y,
   0,0,0,0,0
};

const char dna_nt[128] = {
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,NT_A,0,NT_C,0,0,
   0,NT_G,0,0,0,0,0,0,NT_Z,0,
   0,0,0,0,NT_T,0,0,0,0,0,
   0,0,0,0,0,0,0,NT_A,0,NT_C,
   0,0,0,NT_G,0,0,0,0,0,0,
   NT_Z,0,0,0,0,0,NT_T,0,0,0,
   0,0,0,0,0,0,0,0
};

const char dna_rc[128] = {
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,RC_A,0,RC_C,0,0,
   0,RC_G,0,0,0,0,0,0,NT_Z,0,
   0,0,0,0,RC_T,0,0,0,0,0,
   0,0,0,0,0,0,0,RC_A,0,RC_C,
   0,0,0,RC_G,0,0,0,0,0,0,
   NT_Z,0,0,0,0,0,RC_T,0,0,0,
   0,0,0,0,0,0,0,0
};


// Struct definitions.

typedef struct {
   char * seqname;
   char * chr;
   char * cigar;
   short  flag;
   short  mapq;
   int    locus;
} sam_t;

typedef struct {
   int  pos;
   int  size;
   int  val[];
} stack_t;

typedef struct {
   int    pos;
   int    size;
   char * locus[];
} cstack_t;

typedef struct {
   char   * seqname;
   char   * seq;
   size_t   seqlen;
} ref_t;

typedef struct {
   int     pos;
   int     size;
   ref_t * ref[];
} refstack_t;

// Function headers.
stack_t      * stack_new        (int size);
stack_t      * stack_push       (stack_t ** stackp, int val);
int            stack_pop        (stack_t * stack);
cstack_t     * cstack_new       (int size);
cstack_t     * cstack_push      (cstack_t ** stackp, char * ptr);
char         * cstack_pop       (cstack_t * stack);
refstack_t   * refstack_new     (int size);
refstack_t   * refstack_push    (refstack_t ** stackp, ref_t * ptr);
ref_t        * refstack_pop     (refstack_t * stack);
void           digest_sequence  (ref_t * ref, char * pattern, stack_t ** stack);
refstack_t   * read_genome      (FILE * fg);

// Source.

int main(int argc, char *argv[])
{

   // 1. Parse arguments.
   //    arg list:
   //    1. Organism.
   //    2. RE name.
   //    3. RE sequence.
   //    4. 5' cut offset (from first 5' nucleotide).
   //    5. 3' cut offset (from first 5' nucleotide).

   char * organism;
   char * re_name_;
   char * re_seq;
   int    cut_fw;
   int    cut_rv;

   if (argc > 1 && strcmp(argv[1], "-h") == 0) {
      fprintf(stderr, "%s", HELP_MSG);
      exit(0);
   }
   
   // Parse params.
   if (argc != 6) {
      fprintf(stderr, "usage: %s <organism name> <RE name> <re_sequence> <cut_fw> <cut_rv>\n", argv[0]);
      fprintf(stderr, "type \"%s -h\" for help.\n", argv[0]);
      exit(1);
   }

   organism = argv[1];
   re_name_ = argv[2];
   re_seq = argv[3];
   cut_fw = atoi(argv[4]);
   cut_rv = atoi(argv[5]);

   // Open genome file.
   char * genomepath = malloc(strlen(organism)+17);
   sprintf(genomepath,"db/%s/genome.fasta",organism);

   FILE * genfile = fopen(genomepath,"r");
   if (genfile == NULL) {
      fprintf(stderr, "could not open: %s. Did you create a folder for this organism?\nThe genome of %s must be readable in %s.\n", genomepath, organism, genomepath);
      exit(1);
   }

   // Read genome.
   fprintf(stderr, "reading genome...");
   refstack_t * chr_stack = read_genome(genfile);
   fprintf(stderr, "\tdone\n");
   fclose(genfile);

   // Create new digest file.
   char * re_name = strdup(re_name_);
   lowercase(re_name);
   char * db_path = malloc(strlen(re_name)+strlen(organism)+9);
   sprintf(db_path, "db/%s/%s.isd", organism, re_name);

   struct stat st;
   if (stat(db_path,&st) != -1) {
      fprintf(stderr,"The digestion already exists (File exists: '%s'). Manually remove file to digest again.\n", db_path);
      exit(1);
   }

   int dbfd = open(db_path, O_WRONLY | O_CREAT, 0644);
   if (dbfd < 0) {
      fprintf(stderr, "error while opening: %s.\n",db_path);
      exit(1);
   }
   
   // Write RE seq and cut sites.
   write(dbfd, re_seq, strlen(re_seq)+1);
   write(dbfd, &cut_fw, sizeof(int));
   write(dbfd, &cut_rv, sizeof(int));
   // Write number of chromosomes.
   write(dbfd, &(chr_stack->pos), sizeof(int));
   // Digest genome.
   stack_t * re_sites = stack_new(1024);
   for (int i = 0; i < chr_stack->pos; i++) {
      // Get next chromosome.
      ref_t * ref = chr_stack->ref[i];
      // Find RE sites.
      re_sites->pos = 0;
      fprintf(stderr,"digesting %s...",ref->seqname);
      digest_sequence(ref, re_seq, &re_sites);
      // Write to database. (use low-level write instead of fprint)
      fprintf(stderr,"done\nwrite digestion...");
      // 1. Write chromosome name.
      write(dbfd, ref->seqname, strlen(ref->seqname)+1);
      // 2. Write RE site count.
      write(dbfd, &(re_sites->pos), sizeof(int));
      // 3. Write RE sites.
      size_t offset = 0;
      ssize_t b;
      while ((b = write(dbfd, ((char *)&(re_sites->val))+offset, (re_sites->pos*sizeof(int))-offset)) > 0) offset += b;
      fprintf(stderr,"%ld/%ld bytes written (%d sites)\n",offset,re_sites->pos*sizeof(int), re_sites->pos);
   }

   // Close files and free.
   close(dbfd);
   free(genomepath);
   free(re_name);
   free(re_sites);
   free(db_path);

   return 0;
}


void
digest_sequence
(
 ref_t    * ref,
 char     * pattern,
 stack_t ** stack
)
{
   // Find pattern.
   int pattlen = strlen(pattern);
   stack_push(stack,0);
   for (int i = 0; i < ref->seqlen; i++) {
      int j = 0;
      while (j <= pattlen && (dna_nt[ref->seq[i + j]] & re_nt[pattern[j++]]));
      if (j == pattlen + 1) {
         stack_push(stack, i);
      }
   }
   stack_push(stack,ref->seqlen-1);
}

stack_t *
stack_new
(
 int size
)
{
   stack_t * stack = malloc(sizeof(stack_t) + size*sizeof(int));
   if (!stack)
      return NULL;
   stack->size = size;
   stack->pos = 0;
   return stack;
}

stack_t *
stack_push
(
 stack_t ** stackp,
 int        val
)
{
   stack_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      int newsize = 2*stack->size;
      *stackp = stack = realloc(stack, sizeof(stack_t) + newsize*sizeof(int));
      if (!stack) return NULL;
      stack->size = newsize;
   }

   stack->val[stack->pos++] = val;

   return stack;
}

int
stack_pop
(
 stack_t * stack
)
{
   if (stack->pos == 0) return -1;
   else return stack->val[--stack->pos];
}

cstack_t *
cstack_new
(
 int size
)
{
   cstack_t * stack = malloc(sizeof(cstack_t) + size*sizeof(char *));
   if (!stack)
      return NULL;
   stack->size = size;
   stack->pos = 0;
   return stack;
}

cstack_t *
cstack_push
(
 cstack_t ** stackp,
 char      * ptr
)
{
   cstack_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      int newsize = 2*stack->size;
      *stackp = stack = realloc(stack, sizeof(cstack_t) + newsize*sizeof(char *));
      if (!stack) return NULL;
      stack->size = newsize;
   }

   stack->locus[stack->pos++] = ptr;

   return stack;
}

char *
cstack_pop
(
 cstack_t * stack
)
{
   if (stack->pos == 0) return NULL;
   else return stack->locus[--stack->pos];
}

refstack_t *
refstack_new
(
 int size
)
{
   refstack_t * stack = malloc(sizeof(refstack_t) + size*sizeof(ref_t *));
   if (!stack)
      return NULL;
   stack->size = size;
   stack->pos = 0;
   return stack;
}

refstack_t *
refstack_push
(
 refstack_t ** stackp,
 ref_t       * ptr
)
{
   refstack_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      int newsize = 2*stack->size;
      *stackp = stack = realloc(stack, sizeof(refstack_t) + newsize*sizeof(ref_t *));
      if (!stack) return NULL;
      stack->size = newsize;
   }

   stack->ref[stack->pos++] = ptr;

   return stack;
}

ref_t *
refstack_pop
(
 refstack_t * stack
)
{
   if (stack->pos == 0) return NULL;
   else return stack->ref[--stack->pos];
}

refstack_t *
read_genome
(
 FILE      * fg
)
{
   // Create hash table.
   refstack_t * chrstack = refstack_new(128);
   
   size_t lines = 0, n = 100;
   ssize_t bytes = 0, bufsize;
   char * line = malloc(n);
   ref_t * ref = NULL;
   ENTRY * item;
   while((bytes = getline(&line, &n, fg)) > 0) {
      if (line[bytes-1] == '\n') line[--bytes] = 0;
      // New chromosome.
      if (line[0] == '>') {
         // Realloc current buffer.
         if (ref != NULL) {
            ref->seq = realloc(ref->seq, ref->seqlen*sizeof(char));
         }
         // Crop line.
         strtok(line+1, " ");
         // Create new buffer for the next chromosome.
         bufsize = 1024;
         ref = malloc(sizeof(ref_t));
         ref->seqname = strdup(line+1);
         ref->seq = malloc(bufsize);
         ref->seqlen = 0;
         // Insert buffer into hash table.
         refstack_push(&chrstack, ref);
      }
      // Input is a sequence.
      else {
         // Realloc buffer if full.
         while (ref->seqlen + bytes > bufsize) {
            bufsize *= 2;
            ref->seq = realloc(ref->seq, bufsize*sizeof(char));
            if (ref->seq == NULL) {
               fprintf(stderr, "Error: realloc buffer\n");
               exit(1);
            }
         }
         // Copy data.
         memcpy(ref->seq + ref->seqlen,line,bytes);
         ref->seqlen += bytes;
      }
   }
   free(line);

   if (ref != NULL) {
      ref->seq = realloc(ref->seq, ref->seqlen*sizeof(char));
   }

   
   // Return hash table.
   return chrstack;
}
