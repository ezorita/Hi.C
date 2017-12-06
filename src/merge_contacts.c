#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct {
   char seqname[512];
   int  rev_a;
   char chr_a[512];
   long loc_a;
   int  rev_b;
   char chr_b[512];
   long loc_b;
   int  map_a;
   int  map_b;
} contact_t;


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
   strncpy(cont->seqname, strtok(line, " "), 511);
   cont->rev_a   = strcmp(strtok(NULL, " "), "0") == 0;
   strncpy(cont->chr_a,   strtok(NULL, " "), 511); 
   cont->loc_a   =   atol(strtok(NULL, " "));
                          strtok(NULL, " "); // Dunno what's this field.
   cont->rev_b   = strcmp(strtok(NULL, " "), "0") == 0;
   strncpy(cont->chr_b,   strtok(NULL, " "), 511);
   cont->loc_b   =   atol(strtok(NULL, " "));
                          strtok(NULL, " "); // Still dunno.
   cont->map_a   =   atoi(strtok(NULL, " "));
   cont->map_b   =   atoi(strtok(NULL, " "));
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

   contact_t * cont = calloc(1, sizeof(contact_t));
   contact_t * last = calloc(1, sizeof(contact_t));
   contact_t * tmp;

   int count = 1;

   bytes = getline(&line, &bufsize, fin);

   // First line.
   if (bytes < 1) {
      fprintf(stderr,"error: input file is empty.\n");
      exit(1);
   }
   parse_contact(line, last);

   // Other lines.
   while ((bytes = getline(&line, &bufsize, fin)) > 0) {
      // Parse contact.
      parse_contact(line, cont);
      if (   cont->loc_a == last->loc_a &&
             cont->loc_b == last->loc_b &&
             strcmp(cont->chr_b, last->chr_b) == 0 &&
             strcmp(cont->chr_a, last->chr_a) == 0   ) {
         count++; // Duplicate.
      } else {
         fprintf(stdout, "%s\t%ld\t%s\t%ld\t%d\n",
           last->chr_a, last->loc_a, last->chr_b, last->loc_b, count);
         count = 1;
      }
      tmp  = cont;
      cont = last;
      last = tmp;
   }

   free(cont);
   free(last);
   return 0;

}
