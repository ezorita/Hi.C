#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct {
   char * seqname;
   char * chr_a;
   long   loc_a;
   int    rev_a;
   int    map_a;
   long   beg_a;
   long   end_a;
   char * chr_b;
   long   loc_b;
   int    rev_b;
   int    map_b;
   long   beg_b;
   long   end_b;
} contact_t;

contact_t *
new_contact
(
 void
)
{
   contact_t * cnt = malloc(sizeof(contact_t));
   if (!cnt)
      return NULL;
   cnt->seqname = malloc(500);
   cnt->chr_a   = malloc(500);
   cnt->chr_b   = malloc(500);
   if (!cnt->seqname || !cnt->chr_a || !cnt->chr_b)
      return NULL;
   return cnt;
}

void
free_contact
(
 contact_t * cnt
)
{
   if (cnt == NULL)
      return;
   if (cnt->seqname) {
      free(cnt->seqname);
      cnt->seqname = NULL;
   }
   if (cnt->chr_a) {
      free(cnt->chr_a);
      cnt->chr_a = NULL;
   }
   if (cnt->chr_b) {
      free(cnt->chr_b);
      cnt->chr_b = NULL;
   }
   free(cnt);
}

void
parse_contact
(
 char      * line,
 contact_t * cont
)
{
   size_t len = strlen(line);
   // rstrip line.
   if (line[len-1] == '\n') line[--len] = 0;
   // parse contact.
   strcpy(cont->seqname, strtok(line, "\t"));
   strcpy(cont->chr_a,   strtok(NULL, "\t")); 
   cont->loc_a   = atol(strtok(NULL, "\t"));
   cont->rev_a   = strcmp(strtok(NULL,"\t"), "-") == 0;
   cont->map_a   = atoi(strtok(NULL, "\t"));
   cont->beg_a   = atol(strtok(NULL, "\t"));
   cont->end_a   = atol(strtok(NULL, "\t"));
   strcpy(cont->chr_b,  strtok(NULL, "\t"));
   cont->loc_b   = atol(strtok(NULL, "\t"));
   cont->rev_b   = strcmp(strtok(NULL,"\t"), "-") == 0;
   cont->map_b   = atoi(strtok(NULL, "\t"));
   cont->beg_b   = atol(strtok(NULL, "\t"));
   cont->end_b   = atol(strtok(NULL, "\t"));
}


int main(int argc, char *argv[])
{
   if (argc != 2) {
      fprintf(stderr, "usage: %s file.hcf\n", argv[0]);
      exit(1);
   }
   
   FILE * fin = fopen(argv[1], "r");
   if (!fin) {
      fprintf(stderr, "error opening file: %s\n", argv[1]);
      exit(1);
   }
   
   // Read lines.
   size_t lines = 0;
   size_t bufsize = 200;
   char * line = malloc(bufsize);
   ssize_t bytes = 0;

   contact_t * cont = new_contact();
   contact_t * last = new_contact();
   contact_t * tmp;

   int count = 1;

   bytes = getline(&line, &bufsize, fin);
   // First line..
   if (bytes < 1) {
      fprintf(stderr,"error: input file is empty.\n");
      exit(1);
   }
   parse_contact(line, last);

   // Parse lines.
   while ((bytes = getline(&line, &bufsize, fin)) > 0) {
      // Parse contact.
      parse_contact(line, cont);
      if (cont->beg_b == last->beg_b &&
          cont->beg_a == last->beg_a &&
          strcmp(cont->chr_b, last->chr_b) == 0 &&
          strcmp(cont->chr_a, last->chr_a) == 0
          ) {
         if (cont->loc_a != last->loc_a || cont->loc_b != last->loc_b)
            count++;
      } else {
         fprintf(stdout, "%s\t%ld\t%ld\t%s\t%ld\t%ld\t%d\n",
                 last->chr_a,
                 last->beg_a,
                 last->end_a,
                 last->chr_b,
                 last->beg_b,
                 last->end_b,
                 count
                 );
         count = 1;
      }
      tmp  = cont;
      cont = last;
      last = tmp;
   }

   free_contact(cont);
   free_contact(last);
   return 0;
}


