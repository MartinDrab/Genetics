/* The MIT License

   Copyright (c) 2015 Genome Research Ltd.

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

#include <stdio.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>

int main(int argc, char **argv)
{
    if ( argc<3 )
    {
        fprintf(stderr,
                "Usage:\n"
                "   test-hts-open <output> <mode>\n"
                "Examples:\n"
                "   test-hts-open rmme.bcf wb\n"
                "   test-hts-open rmme.bcf wb0\n"
                "   test-hts-open rmme.bcf wu\n"
                "   test-hts-open rmme.bcf wz\n"
               );
        return 1;
    }

    char *fname = argv[1];
    char *mode  = argv[2];
    bcf_hdr_t *hdr = bcf_hdr_init("w");

    htsFile *fp = hts_open(fname,mode);
    bcf_hdr_write(fp, hdr);
    hts_close(fp);

    return 0;
}

