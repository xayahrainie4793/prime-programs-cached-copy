#ifndef PRIMEGEN_H
#define PRIMEGEN_H

#include <stdint.h>

#define uint32 uint32_t
#define uint64 uint64_t
#define int32 int32_t
#define int64 int64_t

#define PRIMEGEN_WORDS 5000

typedef struct {
  uint32 buf[16][PRIMEGEN_WORDS];
  uint64 p[512]; /* p[num-1] ... p[0], in that order */
  int num;
  int pos; /* next entry to use in buf; WORDS to restart */
  uint64 base;
  uint64 L;
} primegen;

extern void primegen_sieve(primegen *);
extern void primegen_fill(primegen *);

extern void primegen_init(primegen *);
extern uint64 primegen_next(primegen *);
extern uint64 primegen_peek(primegen *);
extern uint64 primegen_count(primegen *,uint64 to);
extern void primegen_skipto(primegen *,uint64 to);

#endif
