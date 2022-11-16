/* files.c -- (C) Geoffrey Reynolds, May 2006.

   Routines for reading and writing sieve and sequence files.

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
#include "sr2sieve.h"

uint32_t read_abcd_file(const char *file_name);
uint32_t read_pfgw_file(const char *file_name);

static int line_counter = 0;
void line_error(const char *msg, const char *fmt, const char *filename)
{
  if (fmt != NULL)
    error("Line %d: %s in %s format file `%s'.",line_counter,msg,fmt,filename);
  else
    error("Line %d: %s in file `%s'.",line_counter,msg,filename);
}

/* Reads a line from file, ignoring blank lines and comments. Returns a
   pointer to the first non-whitespace character in the line.
*/
#define MAX_LINE_LENGTH 256
#define COMMENT_CHAR '#'
const char *read_line(FILE *file)
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

uint32_t read_terms_file(const char *file_name)
{
  FILE *file;
  const char *line;
  int is_abcd = 0, is_pfgw = 0;

  file = xfopen(file_name, "r", error);
  if ((line = read_line(file)) == NULL)
  {
     line_error("Read error", NULL, file_name);
     fclose(file);
     return 0;
  }

  if (!memcmp(line, "ABCD ", 5))
     is_abcd = 1;
  if (!memcmp(line, "ABC ", 4))
     is_pfgw = 1;

  fclose(file);

  if (is_abcd)
     return read_abcd_file(file_name);
  if (is_pfgw)
     return read_pfgw_file(file_name);

  line_error("Unknown terms file format", NULL, file_name);
  return 0;
}

/* This function cannot parse a general ABCD format file, only the specific
   type (ABCD k*b^n+c with fixed k,b,c and variable n) written by srsieve.
*/
uint32_t read_abcd_file(const char *file_name)
{
  FILE *file;
  const char *line;
  uint64_t p;
  uint64_t k;
  uint32_t seq, n, delta, b, n_count, s_count;
  int32_t c;
#if DUAL
  int dual_decided = dual_opt; /* 0 = undecided, 1 = decided */
#endif

  file = xfopen(file_name, "r", error);
  n_count = s_count = seq = n = 0;
  while ((line = read_line(file)) != NULL)
  {
    if (sscanf(line, "%" SCNu32, &delta) == 1)
    {
      if (s_count == 0)
      {
        fclose(file);
        return 0;
      }
      n += delta;
      add_seq_n(seq,n);
      n_count++;
    }
    else
      switch (sscanf(line, "ABCD %" SCNu64 "*%" SCNu32 "^$a%" SCNd32
                     " [%" SCNu32 "] // Sieved to %" SCNu64,
                     &k, &b, &c, &n, &p))
      {
        case 5:
          if (p_min == 0)
            p_min = p;
          else if (p_min != p)
            warning("--pmin=%"PRIu64" from command line overrides pmin=%"
                    PRIu64" from `%s'", p_min, p, file_name);

          /* fall through */

        case 4:
#if (BASE == 0)
          if (b < 2)
            line_error("Invalid base","ABCD",file_name);
          else if (b_term == 0)
            b_term = b;
          else if (b_term != b)
            line_error("Mismatched base","ABCD",file_name);
#else
          if (b != BASE)
            line_error("Invalid base","ABCD",file_name);
#endif
#if DUAL
          if (!dual_decided)
          {
            if (c != 1 && c != -1)
              dual_opt = 1;
            else
              dual_opt = 0;
            dual_decided = 1;
          }
          if (dual_opt)
          {
            if (k != 1)
              line_error("Not dual form b^n+/-k","ABCD",file_name);
            if (c > 0)
              seq = new_seq(c,1);
            else
              seq = new_seq(-c,-1);
          }
          else
#endif
          {
            if (c != 1 && c != -1)
              line_error("Not standard form k*b^n+/-1","ABCD",file_name);
            seq = new_seq(k,c);
          }
          s_count++;
          add_seq_n(seq,n);
          n_count++;
          break;

        default:
          if (s_count == 0)
          {
            fclose(file);
            return 0;
          }
          else
            line_error("Malformed line","ABCD",file_name);
      }
  }

  if (ferror(file))
    line_error("Read error","ABCD",file_name);
  fclose(file);
  report(1,"Read %"PRIu32" term%s for %"PRIu32" sequence%s from ABCD format "
       "file `%s'.",n_count,plural(n_count),s_count,plural(s_count),file_name);

  return n_count;
}

/* This function will read an ABC format file (i.e. .pfgw file). It can handle variable
   or non variable k, b, n, c, but obviously b must be the same throughout. It can also
   parse dual format, with |c|>1, where again k, b, n and c can be variables, but b must
   be the same, and k must be 1.
*/
#define SEQ_ARRAY_GROW_SIZE 64
uint32_t read_pfgw_file(const char *file_name)
{
  FILE *file;
  const char *line;
  uint64_t p;
  uint64_t k, variables[4];
  uint32_t *seqs, *temp_seqs, n, b, n_count, s_count, hash_value, seqs_array_size;
  int32_t c;

  int i, j;

  char *k_str=xmalloc(256),   //To store the appropriate string from the header
       *b_str=xmalloc(256),
       *n_str=xmalloc(256),
       *c_str=xmalloc(256);

  int k_var=0, //Determine the order of variables read in
      b_var=0,
      n_var=0,
      c_var=0,

      p_set=1,    //Flag to check whether p has been read
      num_vars=0, //The number of variables there should be on each line
      num_scn_vars=0,   //The number of variables that are read in on each line
      n_added = 0;   //Flag to check if the current line has been added
#if DUAL
  int dual_decided = dual_opt;
#endif

  for(i=0; i<4; i++)
    variables[i] = 0;
  k = b = n = c = n_count = s_count = 0;
  temp_seqs = seqs = NULL;

  file = xfopen(file_name, "r", error);

  line = read_line(file);
  switch(sscanf(line, "ABC %[$a0-9]*%[$ab0-9]^%[$abc0-9]%[-+$abcd0-9] // Sieved to %" SCNu64"", k_str, b_str, n_str, c_str, &p))
  {
    /* If both n and c are variables, they will be read as one string
      so we sperate them here, and then rescan for p */
    case 3: if(n_str[0]=='$')
            {
              sscanf(n_str, "$%*c%s", c_str);
              sscanf(n_str, "%2[$abc]", n_str);
            }
            else
            {
              sscanf(n_str, "%*" SCNu32 "%s", c_str);
              sscanf(n_str, "%[0-9]", n_str);
            }
            p_set = sscanf(line, "ABC %*s // Sieved to %" SCNu64, &p);
            /* fall through */

    case 5: if(p_set)
            {
              if(p_min == 0)
                p_min = p;
              else if(p_min != p)
                warning("--pmin=%"PRIu64" from command line overrides pmin=%"PRIu64" from `%s'", p_min, p, file_name);
            }
            /* fall through */

    case 4: /* Check which of the values in k*b^n+/-c are variable
               Those that are, make sure we get the correct position
               Those that aren't, get the value from the header*/
            if(k_str[0]=='$')
            {
              k_var = k_str[1]-'a'+1;
              if(k_var>1 || k_var<=0)
                line_error("Invalid header: Variables must be contiguous", "PFGW", file_name);
              num_vars++;
            }
            else
              sscanf(k_str, "%" SCNu64, &k);
           
            if(b_str[0]=='$')
            {
              b_var = b_str[1]-'a'+1;
              if(b_var == k_var)
                line_error("Invalid header: each variable must be unique", "PFGW", file_name);
              if(b_var>2 || b_var<=0 || (k_var+1 != b_var))
                line_error("Invalid header: Variables must be contiguous", "PFGW", file_name);
              num_vars++;
            }
            else
            {
              sscanf(b_str, "%" SCNu32, &b);
#if (BASE==0)
              if(b < 2)
                line_error("Invalid base","PFGW",file_name);
              else
                b_term = b;
#else
              if(b != BASE)
                line_error("Invalid base","PFGW",file_name);
#endif
            }

            if(n_str[0]=='$')
            {
              n_var = n_str[1]-'a'+1;
              if(n_var == k_var || n_var == b_var)
                line_error("Invalid header: each variable must be unique", "PFGW", file_name);
              if(n_var>3 || n_var<=0 || ((k_var+1 != n_var) && (b_var+1 != n_var)))
                line_error("Invalid header: Variables must be contiguous", "PFGW", file_name);
              num_vars++;
            }
            else
              sscanf(n_str, "%" SCNu32, &n);

            if(c_str[0]=='d')
              c_str = "$d";
            if(c_str[0]=='$')
            {
              c_var = c_str[1]-'a'+1;
              if(c_var == k_var || c_var == b_var || c_var == n_var)
                line_error("Invalid header: each variable must be unique", "PFGW", file_name);
              if(c_var>4 || c_var<=0 || ((k_var+1 != c_var) && (b_var+1 != c_var) && (n_var+1 != c_var)))
                line_error("Invalid header: Variables must be contiguous", "PFGW", file_name);
              num_vars++;
            }
            else
            {
              sscanf(c_str, "%" SCNd32, &c);
#if DUAL
              if(!dual_decided)
              {
                if(c!=-1 && c!=1)
                  dual_opt = 1;
                else
                  dual_opt = 0;
                dual_decided = 1;
              }
              if(dual_opt)
              {
              /* If k is not 1, and is not a variable, then we cannot do dual sieve */
                if(k != 1 && k != 0)
                  line_error("Not dual form b^n+/-k","PFGW",file_name);
              }
              else
#endif
              {
                if(c!=-1 && c!=1)
                  line_error("Not standard form k*b^n+/-1","PFGW",file_name);
              }
            }
            break;

    default: line_error("Invalid header", "PFGW", file_name);
  }
  /* End of header reading */

  if(num_vars == 0)
    error("The header must have atleast one variable");

  /* If k and c are fixed, then there is only one sequence */
  if(k_var==0 && c_var==0)
  {
    seqs = xmalloc(sizeof(uint32_t));
#if DUAL
    if(dual_opt)
    {
      if(c>0)
        seqs[0] = new_seq(c, 1);
      else
        seqs[0] = new_seq(-c, -1);
    }
    else
#endif
    {
      seqs[0] = new_seq(k, c);
    }
    s_count=1;
  /* Otherwise create the sequence array, with SEQ_ARRAY_GROW_SIZE elements */
  }
  else
  {
    seqs = xmalloc(sizeof(uint32_t)*SEQ_ARRAY_GROW_SIZE);
    temp_seqs = xmalloc(sizeof(uint32_t)*SEQ_ARRAY_GROW_SIZE);
    seqs_array_size = SEQ_ARRAY_GROW_SIZE;
    for(i=0; i<seqs_array_size; i++)
    {
      seqs[i] = 0;
      temp_seqs[i] = 0;
    }
  }

  /* Begin reading of sieve */
  while((line = read_line(file)) != NULL)
  {
    n_added = 0;
    if((num_scn_vars = sscanf(line, "%" SCNu64 " %" SCNu64 " %" SCNu64 " %" SCNu64, &variables[0], &variables[1], &variables[2], &variables[3])) != num_vars)
      line_error("Malformed line", "PFGW", file_name);
    else
    {  /* For all of the non static variables, copy the appropriate value
          that we just read */
      if(k_var!=0)
         k = variables[k_var-1];
#if DUAL
      if(dual_opt && k!=1)
        line_error("Not dual form b^n+/-k","PFGW",file_name);
#endif

      if(b_var!=0)
      {
        b = variables[b_var-1];
#if (BASE==0)
        if(b<2)
          line_error("Invalid base", "PFGW", file_name);
        if(b_term != 0 && b != b_term)
          line_error("Mismatched base", "PFGW", file_name);
        else
          b_term = b;
#else
        if(b != BASE)
          line_error("Invalid base", "PFGW", file_name);
#endif
      }

      if(n_var!=0) {
        if(variables[n_var-1] < n)
          line_error("File is not sorted by n ascending", "PFGW", file_name);
        n = variables[n_var-1];
      }

      if(c_var!=0)
      {
        c = (int32_t) variables[c_var-1];
#if DUAL
        if(!dual_decided)
        {
          if(c!=-1 && c!=1)
            dual_opt = 1;
          else
            dual_opt = 0;
          dual_decided = 1;
        }
        if(dual_opt)
        {
          if(k!=1)
            line_error("Not dual form b^n+/-k","PFGW",file_name);
        }
        else
#endif
        {
          if(c!=-1 && c!=1)
            line_error("Not standard form k*b^n+/-1","PFGW",file_name);
        }
      }

      /* If k and c are static, add the current value to the sequence */
      if(k_var==0 && c_var==0)
      {
        add_seq_n(seqs[0], n);
        n_count++;

      /* Otherwise, create the hash, and search the sequences array
         for the current sequence. When it is found, add the current
         value to it. If it isn't found, then create it, and add the
         current value */
      }
      else
      {
#if DUAL
        if(dual_opt)
        {
          if(c>0)
            hash_value = c%seqs_array_size;
          else
            hash_value = (-c)%seqs_array_size;
        }
        else
#endif
        {
          hash_value = (k+c)%seqs_array_size;
        }

        for(i=hash_value; i<seqs_array_size; i++)
        {
          if(seqs[i]==0)
          {
#if DUAL
            if(dual_opt)
            {
              if(c>0)
                seqs[i] = new_seq(c, 1)+1;
              else
                seqs[i] = new_seq(-c, -1)+1;
            }
       else
#endif
            {
              seqs[i] = new_seq(k, c)+1;
            }
            s_count++;
            add_seq_n(seqs[i]-1, n);
            n_count++;
            n_added = 1;
            break;
          }
#if DUAL
          else if((dual_opt) && ((SEQ[seqs[i]-1].k == c && SEQ[seqs[i]-1].c == 1) || (SEQ[seqs[i]-1].k == -c && SEQ[seqs[i]-1].c == -1)))
          {
            if(SEQ[seqs[i]-1].N[SEQ[seqs[i]-1].ncount-1] == n)
              line_error("Repeated n for 1 sequence", "PFGW", file_name);
            add_seq_n(seqs[i]-1, n);
            n_count++;
            n_added = 1;
            break;
          }
#endif
          else if(SEQ[seqs[i]-1].k == k && SEQ[seqs[i]-1].c == c)
          {
            if(SEQ[seqs[i]-1].N[SEQ[seqs[i]-1].ncount-1] == n)
              line_error("Repeated n for 1 sequence", "PFGW", file_name);
            add_seq_n(seqs[i]-1, n);
            n_count++;
            n_added = 1;
            break;
     }
        } /* End of first for loop */

        /* If we haven't found the sequence, or a space, by the end
           of the array, start at the beginning, and continue upwards
           to the hash value */
        if(n_added == 0)
          for(i=0; i<hash_value; i++)
          {
            if(seqs[i]==0)
            {
#if DUAL
              if(dual_opt)
              {
                if(c>0)
                  seqs[i] = new_seq(c, 1)+1;
                else
                  seqs[i] = new_seq(-c, -1)+1;
              }
              else
#endif
              {
                seqs[i] = new_seq(k,c)+1;
              }
              s_count++;
              add_seq_n(seqs[i]-1, n);
              n_count++;
              n_added = 1;
              break;
            }
#if DUAL
            else if((dual_opt) && ((SEQ[seqs[i]-1].k == c && SEQ[seqs[i]-1].c == 1) || (SEQ[seqs[i]-1].k == -c && SEQ[seqs[i]-1].c == -1)))
            {
              if(SEQ[seqs[i]-1].N[SEQ[seqs[i]-1].ncount-1] == n)
                line_error("Repeated n for 1 sequence", "PFGW", file_name);
              add_seq_n(seqs[i]-1, n);
              n_count++;
              n_added = 1;
              break;
            }
#endif
            else if(SEQ[seqs[i]-1].k == k && SEQ[seqs[i]-1].c == c)
            {
              if(SEQ[seqs[i]-1].N[SEQ[seqs[i]-1].ncount-1] == n)
                line_error("Repeated n for 1 sequence", "PFGW", file_name);
              add_seq_n(seqs[i]-1, n);
              n_count++;
              n_added = 1;
              break;
            }
          }
      } /* End of adding current line */
    } /* End of processing current line */

    /* If the hash table's (sequences array) load factor increases above 0.7,
       then we want to increase its size. Create a temp copy of the sequences,
       Increase the sequence array size by SEQ_ARRAY_GROW_SIZE. Then rehash
       each of the non 0 values on the temp array, and store them in the array.
       finally, resize the temp array, for the next reallocation. */
    if(s_count > seqs_array_size*0.7 && (k_var!=0 || c_var!=0))
    {
      for(i=0; i<seqs_array_size; i++)
        temp_seqs[i] = seqs[i];
      seqs_array_size+=SEQ_ARRAY_GROW_SIZE;
      seqs = xrealloc(seqs, seqs_array_size*sizeof(uint32_t));
      for(i=0; i<seqs_array_size; i++)
        seqs[i]=0;

      for(i=0; i<seqs_array_size-SEQ_ARRAY_GROW_SIZE; i++)
        if(temp_seqs[i]!=0)
        {
          n_added = 0;
#if DUAL
          if(dual_opt)
          {
            hash_value = (SEQ[temp_seqs[i]-1].k)%seqs_array_size;
          }
          else
#endif
          {
            hash_value = (SEQ[temp_seqs[i]-1].k + SEQ[temp_seqs[i]-1].c)%seqs_array_size;
          }
          for(j=hash_value; j<seqs_array_size; j++)
            if(seqs[j]==0)
            {
              seqs[j] = temp_seqs[i];
              n_added = 1;
              break;
            }

         if(n_added == 0)
           for(j=0; j<hash_value; j++)
             if(seqs[j]==0)
             {
               seqs[j] = temp_seqs[i];
               n_added = 1;
               break;
             }
        }

        temp_seqs = xrealloc(temp_seqs, seqs_array_size*sizeof(uint32_t));
        for(i=0; i<seqs_array_size; i++)
          temp_seqs[i] = 0;
    } /* End of hash table reallocation */
  } /* End of file processing */

  if (ferror(file))
    line_error("Read error","PFGW",file_name);
  fclose(file);
  report(1,"Read %"PRIu32" term%s for %"PRIu32" sequence%s from PFGW format file `%s'.",n_count,plural(n_count),s_count,plural(s_count),file_name);

  return n_count;
}

#if (BASE==5)
void overwrite_abcd(uint64_t p, const char *file_name, uint64_t del_k)
{
  FILE *file;
  uint32_t i, seq, n, s_count, t_count;

  if ((file = xfopen(file_name, "w", warning)) == NULL)
    return;

  s_count = t_count = 0;
  for (seq = 0; seq < seq_count; seq++)
    if (SEQ[seq].k == del_k)
    {
      report(1,"Removed all %"PRIu32" terms for sequence %s from the sieve.",
             SEQ[seq].ncount, seq_str(seq));
    }
    else
    {
      fprintf(file, "ABCD %" PRIu64 "*%u^$a%+" PRId32,
              SEQ[seq].k, (unsigned int)BASE, SEQ[seq].c);
      n = SEQ[seq].N[0];
      if (p > 0)
      {
        fprintf(file," [%" PRIu32 "] // Sieved to %" PRIu64 " with srsieve\n",
                n,p);
        p = 0;
      }
      else
        fprintf(file," [%" PRIu32 "]\n", n);

      for (i = 1; i < SEQ[seq].ncount; i++)
      {
        fprintf(file,"%" PRIu32 "\n",SEQ[seq].N[i]-n);
        n = SEQ[seq].N[i];
      }
      t_count += i;
      s_count ++;
    }

  report(1,"Wrote %"PRIu32" term%s for %"PRIu32" sequence%s to ABCD format "
       "file `%s'.",t_count,plural(t_count),s_count,plural(s_count),file_name);

  xfclose(file,file_name);
}

void write_newpgen_k(uint32_t seq, uint32_t n0, uint32_t n1)
{
  FILE *file;
  char file_name[FILENAME_MAX+1];
  int i, count;

  assert(ABS(SEQ[seq].c) == 1);

  snprintf(file_name, FILENAME_MAX, "%"PRIu64".txt", SEQ[seq].k);
  file_name[FILENAME_MAX] = '\0';
  file = xfopen(file_name, "w", error);
  if (SEQ[seq].c == 1)
    fprintf(file, "%"PRIu64":P:1:%u:257\n", p_min, (unsigned int)BASE);
  else
    fprintf(file, "%"PRIu64":M:1:%u:258\n", p_min, (unsigned int)BASE);

  for (i = 0, count = 0; i < SEQ[seq].ncount && SEQ[seq].N[i] < n1; i++)
    if (SEQ[seq].N[i] >= n0)
    {
      fprintf(file, "%"PRIu64" %"PRIu32"\n", SEQ[seq].k, SEQ[seq].N[i]);
      count++;
    }

  report(1,"Wrote %d term%s for sequence %s to %s format file `%s'.",
         count, plural(count), seq_str(seq), "NewPGen", file_name);
  xfclose(file,file_name);
}

void write_abc_k(uint32_t seq, uint32_t n0, uint32_t n1)
{
  FILE *file;
  char file_name[FILENAME_MAX+1];
  int i, count;

  assert(ABS(SEQ[seq].c) == 1);

  snprintf(file_name, FILENAME_MAX, "%"PRIu64".abc", SEQ[seq].k);
  file_name[FILENAME_MAX] = '\0';
  file = xfopen(file_name, "w", error);
  fprintf(file, "ABC %"PRIu64"*%u^$a%+"PRId32"\n",
          SEQ[seq].k, (unsigned int)BASE, SEQ[seq].c);
  for (i = 0, count = 0; i < SEQ[seq].ncount && SEQ[seq].N[i] < n1; i++)
    if (SEQ[seq].N[i] >= n0)
    {
      fprintf(file, "%"PRIu32"\n", SEQ[seq].N[i]);
      count++;
    }

  report(1,"Wrote %d term%s for sequence %s to %s format file `%s'.",
         count, plural(count), seq_str(seq), "ABC", file_name);
  xfclose(file,file_name);
}
#endif /* BASE==5 */

uint64_t read_checkpoint(uint64_t pmin, uint64_t pmax)
{
  FILE *file;
  uint64_t p;
  double cpu_secs = 0.0, elapsed_secs = 0.0;
  uint32_t count = 0;

  if ((file = fopen(checkpoint_file_name,"r")) == NULL)
    return pmin;

  if (fscanf(file,"pmin=%"SCNu64",factor_count=%"SCNu32
             ",cpu_secs=%lf,frac_done=%*f,elapsed_secs=%lf",
             &p,&count,&cpu_secs,&elapsed_secs) < 1)
    error("Cannot read checkpoint from `%s'.",checkpoint_file_name);

  xfclose(file,checkpoint_file_name);

  if (p > pmin && p < pmax)
  {
#if SOBISTRATOR_OPT
    if (!sobistrator_opt)
#endif
      report(1,"Resuming from checkpoint pmin=%" PRIu64 " in `%s'.",
             p, checkpoint_file_name);

    factor_count = count;
    set_accumulated_cpu(cpu_secs);
    set_accumulated_time(elapsed_secs);
    return p;
  }

  return pmin;
}

void write_checkpoint(uint64_t p)
{
  FILE *file;

  if ((file = xfopen(checkpoint_file_name,"w",warning)) != NULL)
  {
    fprintf(file,"pmin=%"PRIu64",factor_count=%"PRIu32
            ",cpu_secs=%.3f,frac_done=%f,elapsed_secs=%.3f\n",
            p,factor_count,get_accumulated_cpu(),frac_done(p),
            get_accumulated_time());
    xfclose(file,checkpoint_file_name);
  }
}

void remove_checkpoint(void)
{
  remove(checkpoint_file_name);
}

#if (BASE==2 || BASE==0)
uint32_t read_dat_file(const char *file_name, int32_t c)
{
  FILE *file;
  const char *line;
  uint32_t seq, k, n, delta, n_count, s_count;

#if (BASE==0)
  assert(b_term == 0 || b_term == 2);
  b_term = 2;
#endif

  file = xfopen(file_name,"r",error);
  n_count = s_count = seq = n = 0;
  while ((line = read_line(file)) != NULL)
  {
    if (sscanf(line,"%"SCNu32,&delta) == 1)
    {
      if (s_count > 0)
      {
        if (line[0] == '+')
          n += delta;
        else
          n = delta;
        add_seq_n(seq,n);
        n_count++;
      }
      /* else continue; skip first few unused lines. */
    }
    else if (sscanf(line,"k=%"SCNu32,&k) == 1)
    {
      seq = new_seq(k,c);
      s_count++;
      n = 0;
    }
    else
      line_error("Malformed line","dat",file_name);
  }

  if (ferror(file))
    line_error("Read error","dat",file_name);
  fclose(file);
  printf("Read %"PRIu32" term%s for %"PRIu32" sequence%s from dat format file"
         " `%s'.\n",n_count,plural(n_count),s_count,plural(s_count),file_name);

  return n_count;
}
#endif

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
