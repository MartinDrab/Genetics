/*  test/test-regidx.c -- Regions index test harness.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include <htslib/regidx.h>

static int verbose = 0;

void debug(const char *format, ...)
{
    if ( verbose<2 ) return;
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}
void info(const char *format, ...)
{
    if ( verbose<1 ) return;
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}
void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

int custom_parse(const char *line, char **chr_beg, char **chr_end, reg_t *reg, void *payload, void *usr)
{
    // Use the standard parser for CHROM,FROM,TO
    int i, ret = regidx_parse_tab(line,chr_beg,chr_end,reg,NULL,NULL);
    if ( ret!=0 ) return ret;

    // Skip the fields that were parsed above
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    for (i=0; i<3; i++)
    {
        while ( *ss && !isspace(*ss) ) ss++;
        if ( !*ss ) return -2;  // wrong number of fields
        while ( *ss && isspace(*ss) ) ss++;
    }
    if ( !*ss ) return -2;

    // Parse the payload
    char *se = ss;
    while ( *se && !isspace(*se) ) se++;
    char **dat = (char**) payload;
    *dat = (char*) malloc(se-ss+1);
    memcpy(*dat,ss,se-ss+1);
    (*dat)[se-ss] = 0;
    return 0;
}
void custom_free(void *payload)
{
    char **dat = (char**)payload;
    free(*dat);
}

void test_custom_payload(void)
{
    // Init index with no file name, we will insert the regions manually
    regidx_t *idx = regidx_init(NULL,custom_parse,custom_free,sizeof(char*),NULL);
    if ( !idx ) error("init failed\n");

    // Insert regions
    char *line;
    line = "1 10000000 10000000 1:10000000-10000000"; if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
    line = "1 20000000 20000001 1:20000000-20000001"; if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
    line = "1 20000002 20000002 1:20000002-20000002"; if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
    line = "1 30000000 30000000 1:30000000-30000000"; if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);

    // Finish initialization
    regidx_insert(idx,NULL);

    // Test 
    regitr_t itr;
    int from, to;

    from = to = 10000000;
    if ( !regidx_overlap(idx,"1",from-1,to-1,&itr) ) error("query failed: 1:%d-%d\n",from,to);
    if ( strcmp("1:10000000-10000000",REGITR_PAYLOAD(itr,char*)) ) error("query failed: 1:%d-%d vs %s\n", from,to,REGITR_PAYLOAD(itr,char*));
    if ( !regidx_overlap(idx,"1",from-2,to-1,&itr) ) error("query failed: 1:%d-%d\n",from-1,to);
    if ( !regidx_overlap(idx,"1",from-2,to+3,&itr) ) error("query failed: 1:%d-%d\n",from-1,to+2);
    if ( regidx_overlap(idx,"1",from-2,to-2,&itr) ) error("query failed: 1:%d-%d\n",from-1,to-1);

    from = to = 20000000;
    if ( !regidx_overlap(idx,"1",from-1,to-1,&itr) ) error("query failed: 1:%d-%d\n",from,to);

    from = to = 20000002;
    if ( !regidx_overlap(idx,"1",from-1,to-1,&itr) ) error("query failed: 1:%d-%d\n",from,to);

    from = to = 30000000;
    if ( !regidx_overlap(idx,"1",from-1,to-1,&itr) ) error("query failed: 1:%d-%d\n",from,to);

    // Clean up
    regidx_destroy(idx);
}

void create_line_bed(char *line, char *chr, int start, int end)
{
    sprintf(line,"%s\t%d\t%d\n",chr,start-1,end);
}
void create_line_tab(char *line, char *chr, int start, int end)
{
    sprintf(line,"%s\t%d\t%d\n",chr,start,end);
}
void create_line_reg(char *line, char *chr, int start, int end)
{
    sprintf(line,"%s:%d-%d\n",chr,start,end);
}

typedef void (*set_line_f)(char *line, char *chr, int start, int end);

void test(set_line_f set_line, regidx_parse_f parse)
{
    regidx_t *idx = regidx_init(NULL,parse,NULL,0,NULL);
    if ( !idx ) error("init failed\n");

    char line[250], *chr = "1";
    int i, n = 10, start, end, nhit;
    for (i=1; i<n; i++)
    {
        start = end = 10*i;
        set_line(line,chr,start,end);
        if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
        debug("insert: %s", line);

        start = end = 10*i + 1;
        set_line(line,chr,start,end);
        if ( regidx_insert(idx,line)!=0 ) error("insert failed: %s\n", line);
        debug("insert: %s", line);
    }
    regidx_insert(idx,NULL);

    regitr_t itr;
    for (i=1; i<n; i++)
    {
        // no hit
        start = end = 10*i - 1;
        if ( regidx_overlap(idx,chr,start-1,end-1,&itr) ) error("query failed, there should be no hit: %s:%d-%d\n",chr,start,end);


        // one hit
        start = end = 10*i;
        if ( !regidx_overlap(idx,chr,start-1,end-1,&itr) ) error("query failed, there should be a hit: %s:%d-%d\n",chr,start,end);
        nhit = 0;
        while ( REGITR_OVERLAP(itr,start-1,end-1) )
        {
            nhit++;
            itr.i++;
        }
        if ( nhit!=1 ) error("query failed, expected one hit, found %d: %s:%d-%d\n",nhit,chr,start,end);


        // one hit
        start = end = 10*i+1;
        if ( !regidx_overlap(idx,chr,start-1,end-1,&itr) ) error("query failed, there should be a hit: %s:%d-%d\n",chr,start,end);
        nhit = 0;
        while ( REGITR_OVERLAP(itr,start-1,end-1) )
        {
            nhit++;
            itr.i++;
        }
        if ( nhit!=1 ) error("query failed, expected one hit, found %d: %s:%d-%d\n",nhit,chr,start,end);


        // two hits
        start = 10*i; end = start+1;
        if ( !regidx_overlap(idx,chr,start-1,end-1,&itr) ) error("query failed, there should be a hit: %s:%d-%d\n",chr,start,end);
        nhit = 0;
        while ( REGITR_OVERLAP(itr,start-1,end-1) )
        {
            nhit++;
            itr.i++;
        }
        if ( nhit!=2 ) error("query failed, expected two hits, found %d: %s:%d-%d\n",nhit,chr,start,end);

    }

    regidx_destroy(idx);
}

static void usage(void)
{
    fprintf(stderr, "Usage: test-regidx [OPTIONS]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -h, --help          this help message\n");
    fprintf(stderr, "   -v, --verbose       increase verbosity by giving multiple times\n");

    exit(1);
}

int main(int argc, char **argv)
{
    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"verbose",0,0,'v'},
        {0,0,0,0}
    };
    char c;
    while ((c = getopt_long(argc, argv, "hv",loptions,NULL)) >= 0) 
    {
        switch (c)
        {
            case 'v': verbose++; break;
            default: usage(); break;
        }
    }

    info("Testing TAB\n");
    test(create_line_tab,regidx_parse_tab);

    info("Testing REG\n");
    test(create_line_reg,regidx_parse_reg);

    info("Testing BED\n");
    test(create_line_bed,regidx_parse_bed);

    info("Testing custom payload\n");
    test_custom_payload();

    return 0;
}


