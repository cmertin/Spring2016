
/**
  @file dendro.h
  @author Hari Sundar   hari@cs.utah.edu
  */

#ifndef __DENDRO_H__
#define __DENDRO_H__

// the configured options and settings for Tutorial
#define DENDRO_VERSION_MAJOR 4
#define DENDRO_VERSION_MINOR 0


#ifdef USE_64BIT_INDICES
#define DendroIntL long long
#define DendroIntLSpecifier %lld
#define DendroUIntLSpecifier %llu
#else
#define DendroIntL int
#define DendroIntLSpecifier %d
#define DendroUIntLSpecifier %u
#endif

#endif

