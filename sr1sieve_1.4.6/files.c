/* files.c -- (C) Geoffrey Reynolds, May 2006.

   Routines for reading and writing sieve files.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <assert.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include "sr1sieve.h"

format_t outputFileFormat = FF_ABCD;
uint32_t last_n_term;

void read_abc_file(const char *file_name);
void read_abcd_file(const char *file_name);
void read_newpgen_file(const char *file_name);

void write_abc_file(int scroll, uint64_t p, const char *file_name);
void write_abcd_file(int scroll, uint64_t p, const char *file_name);
void write_newpgen_file(int scroll, uint64_t p, const char *file_name);

static int line_counter = 0;
static void line_error(const char *msg, const char *file_name)
{
  error("Line %d: %s in input file `%s'.",line_counter,msg,file_name);
}

/* Reads a line from file, ignoring blank lines and comments. Returns a
   pointer to the first non-whitespace character in the line.
*/
#define MAX_LINE_LENGTH 256
#define COMMENT_CHAR '#'
static const char *read_line(FILE *file)
{
  static char buf[MAX_LINE_LENGTH];
  const char *ptr;

  while ((ptr = fgets(buf,MAX_LINE_LENGTH,file)) != NULL)
  {
    line_counter++;
    while (isspace(*ptr))
      ptr++;
    if (*ptr != COMMENT_CHAR && *ptr != '\0')
      break;
  }

  return ptr;
}

FILE *xfopen(const char *fn, const char *mode, void (*fun)(const char *,...))
{
  FILE *ret;

  ret = fopen(fn,mode);
  if (ret == NULL && fun != NULL)
    fun("Failed to open %sput file `%s'.", (*mode == 'r') ? "in" : "out", fn);

  line_counter = 0;
  return ret;
}

void xfclose(FILE *file, const char *fn)
{
  if (file != NULL && fclose(file))
    warning("Problem closing file `%s'.", fn);
}

void read_input_file(const char *file_name)
{
  FILE *file;
  const char *line;
  
  file = xfopen(file_name,"r",error);
  line = read_line(file);
  
  if (line == 0)
    line_error("Invalid header",file_name);
 
  fclose(file);
  
  if (!memcmp(line, "ABCD ", 5))
     read_abcd_file(file_name);
  else if (!memcmp(line, "ABC ", 4))
     read_abc_file(file_name);
  else
     read_newpgen_file(file_name);
}

void read_abc_file(const char *file_name)
{
  FILE *file;
  const char *line;
  uint64_t k;
  uint32_t n;
  int32_t  c;

  /* Try to read header */
  file = xfopen(file_name,"r",error);
  line = read_line(file);

  if (sscanf(line, "ABC $a*%" SCNu32"^$b$c // Sieved to %" SCNu64"", &b_term, &k) != 2)
    line_error("Invalid header",file_name);

  if (b_term < 2)
    line_error("Invalid base",file_name);
  if (p_min == 0)
    p_min = k;
  else if (p_min != k)
    warning("--pmin=%"PRIu64" from command line overrides pmin=%"PRIu64
            " from `%s'", p_min, k, file_name);

  if ((line = read_line(file)) == NULL)
    error("Empty sieve file `%s'.",file_name);

  if (sscanf(line, "%" SCNu64 " %" SCNu32 " %" SCNd32 "", &k_term, &n, &c_term) != 3)
    line_error("Malformed line",file_name);
  if (k_term < 2)
    line_error("k term too small",file_name);
  if (c_term == 0)
    line_error("c term cannot be 0",file_name);
  
  add_seq_n(n);

  while ((line = read_line(file)) != NULL)
  {
    if (sscanf(line, "%" SCNu64 " %" SCNu32 " %" SCNd32 "", &k, &n, &c) != 3)
      line_error("Malformed line",file_name);
    if (k != k_term)
      line_error("Invalid k term",file_name);
    if (c != c_term)
      line_error("Invalid c term",file_name);
    if (n <= N[ncount-1])
      line_error("Non-increasing n term",file_name);
    add_seq_n(n);
  }

  if (ferror(file))
    line_error("Read error",file_name);
  fclose(file);
  report(1,"Read %"PRIu32" term%s for %s from ABC file `%s'.",
         ncount, plural(ncount), kbc_str(), file_name);
}

void read_abcd_file(const char *file_name)
{
  FILE *file;
  const char *line;
  uint64_t k;
  uint32_t n, last_n_term;

  /* Try to read header */
  file = xfopen(file_name,"r",error);
  line = read_line(file);
  
  if (sscanf(line,  "ABCD %" SCNu64"*%" SCNu32"^$a%" SCNd32" [%" SCNu32"] // Sieved to %" SCNu64"", &k_term, &b_term, &c_term, &n, &k) != 5)
    line_error("Invalid header",file_name);

  if (b_term < 2)
    line_error("Invalid base",file_name);
  if (p_min == 0)
    p_min = k;
  else if (p_min != k)
    warning("--pmin=%"PRIu64" from command line overrides pmin=%"PRIu64
            " from `%s'", p_min, k, file_name);

  if (k_term < 2)
    line_error("k term too small",file_name);
  if (c_term == 0)
    line_error("c term cannot be 0",file_name);
  
  add_seq_n(n);
  last_n_term = n;

  while ((line = read_line(file)) != NULL)
  {
    if (!memcmp(line, "ABCD ", 5))
      line_error("Cannot support multiple sequences", file_name);
    if (sscanf(line, "%" SCNu32 "", &n) != 1)
      line_error("Malformed line",file_name);
    n += last_n_term;
    if (n <= N[ncount-1])
      line_error("Non-increasing n term",file_name);
    add_seq_n(n);
    last_n_term = n;
  }

  if (ferror(file))
    line_error("Read error",file_name);
  fclose(file);
  report(1,"Read %"PRIu32" term%s for %s from ABCD file `%s'.",
         ncount, plural(ncount), kbc_str(), file_name);
}

void read_newpgen_file(const char *file_name)
{
  FILE *file;
  const char *line;
  char mode_char;
  int mode_bits;
  uint64_t k;
  uint32_t n;

  /* Try to read header */
  file = xfopen(file_name,"r",error);
  line = read_line(file);
  /* headers are of the form <p_start>:<mode_char>:<int>:<base>:<mode_bits>
   */
  if (sscanf(line, "%" SCNu64 ":%c:%*d:%" SCNu32 ":%d", &k, &mode_char, &b_term, &mode_bits) != 4)
    line_error("Invalid header",file_name);

  if (b_term < 2)
    line_error("Invalid base",file_name);
  if (mode_char == 'M' && mode_bits == 258)
    c_term = -1;
  else if (mode_char == 'P' && mode_bits == 257)
    c_term = 1;
  else
    line_error("Wrong sieve type",file_name);
  if (p_min == 0)
    p_min = k;
  else if (p_min != k)
    warning("--pmin=%"PRIu64" from command line overrides pmin=%"PRIu64
            " from `%s'", p_min, k, file_name);

  if ((line = read_line(file)) == NULL)
    error("Empty sieve file `%s'.",file_name);

  if (sscanf(line, "%" SCNu64 " %" SCNu32, &k, &n) != 2)
    line_error("Malformed line",file_name);
  if (k < 2)
    line_error("k term too small",file_name);
  k_term = k;
  add_seq_n(n);

  while ((line = read_line(file)) != NULL)
  {
    if (sscanf(line, "%"SCNu64" %"SCNu32, &k, &n) != 2)
      line_error("Malformed line",file_name);
    if (k != k_term)
      line_error("Invalid k term",file_name);
    if (n <= N[ncount-1])
      line_error("Non-increasing n term",file_name);
    add_seq_n(n);
  }

  if (ferror(file))
    line_error("Read error",file_name);
  fclose(file);
  report(1,"Read %"PRIu32" term%s for %s from NewPGen file `%s'.",
         ncount, plural(ncount), kbc_str(), file_name);
}

static void write_abc_term(uint32_t n, void *file)
{
  fprintf((FILE *)file,"%"PRIu64" %"PRIu32" %"PRId32"\n",k_term,n,c_term);
}

static void write_abcd_term(uint32_t n, void *file)
{
  if (n == last_n_term)
    return;
  
  fprintf((FILE *)file, "%"PRIu32"\n", n - last_n_term);
  last_n_term = n;
}

static void write_newpgen_term(uint32_t n, void *file)
{
  fprintf((FILE *)file,"%"PRIu64" %"PRIu32"\n",k_term,n);
}

void write_output_file(int scroll, uint64_t p, const char *file_name)
{
   if (outputFileFormat == FF_NEWPGEN)
      write_newpgen_file(scroll, p, file_name);
   if (outputFileFormat == FF_ABC)
      write_abc_file(scroll, p, file_name);
   if (outputFileFormat == FF_ABCD)
      write_abcd_file(scroll, p, file_name);
}

void write_abc_file(int scroll, uint64_t p, const char *file_name)
{
  if (file_name != NULL)
  {
    FILE *file;

    if ((file = xfopen(file_name,"w",warning)) != NULL)
    {
      uint32_t count;

      fprintf(file, "ABC $a*%"PRIu32"^$b$c // Sieved to %"PRIu64"\n", b_term, p);

      count = for_each_term(write_abc_term, file);

      xfclose(file,file_name);
      report(scroll,"Wrote %"PRIu32" term%s for %s to ABC file `%s'.",
             count, plural(count), kbc_str(), file_name);
    }
  }
}

void write_abcd_file(int scroll, uint64_t p, const char *file_name)
{
  if (file_name != NULL)
  {
    last_n_term = get_first_term();

    if (last_n_term == 0)
    {
      report(scroll, "Nothing to write as no terms remain.");
      return;
    } 
     
    FILE *file;

    if ((file = xfopen(file_name,"w",warning)) != NULL)
    {
      uint32_t count;

      fprintf(file, "ABCD %"PRIu64"*%"PRIu32"^$a%+d [%"PRIu32"] // Sieved to %"PRIu64"\n", k_term, b_term, c_term, last_n_term, p);

      count = for_each_term(write_abcd_term, file);

      xfclose(file,file_name);
      report(scroll,"Wrote %"PRIu32" term%s for %s to ABCD file `%s'.",
             count, plural(count), kbc_str(), file_name);
    }
  }
}

void write_newpgen_file(int scroll, uint64_t p, const char *file_name)
{
  if (file_name != NULL)
  {
    FILE *file;

    if ((file = xfopen(file_name,"w",warning)) != NULL)
    {
      uint32_t count;

      if (c_term == -1)
        fprintf(file,"%"PRIu64":M:1:%"PRIu32":258\n",p,b_term);
      else
        fprintf(file,"%"PRIu64":P:1:%"PRIu32":257\n",p,b_term);

      count = for_each_term(write_newpgen_term,file);

      xfclose(file,file_name);
      report(scroll,"Wrote %"PRIu32" term%s for %s to NewPGen file `%s'.",
             count, plural(count), kbc_str(), file_name);
    }
  }
}

#if USE_COMMAND_LINE_FILE
#include <string.h>
#include <stdlib.h>
void read_argc_argv(int *argc, char ***argv, const char *file_name)
{
  FILE *file;
  char *line;
  char **strv;
  int i, len;

  if ((file = xfopen(file_name,"r",NULL)) == NULL)
    return;

  line = xmalloc(1024);
  if (fgets(line,1023,file) == NULL)
  {
    free(line);
    fclose(file);
    return;
  }

  len = 16;
  strv = xmalloc(len*sizeof(char *));
  if ((strv[0] = strtok(line," \t\r\n")) == NULL)
  {
    free(strv);
    free(line);
    fclose(file);
    return;
  }

  for (i = 1; ((strv[i] = strtok(NULL," \t\r\n")) != NULL); i++)
    if (len <= i+1)
    {
      len += 16;
      strv = xrealloc(strv,len*sizeof(char *));
    }

  report(1,"(Read %d command line arguments from file `%s').",i-1,file_name);

  *argc = i;
  *argv = strv;
  fclose(file);
}
#endif
