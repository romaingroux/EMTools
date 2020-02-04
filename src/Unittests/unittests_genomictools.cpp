#include <UnitTest++/UnitTest++.h>
#include <seqan/bam_io.h>
#include <GenomeRegion.hpp>
#include <CorrelationMatrixCreator.hpp>
#include <Matrix2D.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <stdexcept>        // std::invalid_argument

std::string file_bed = "/local/groux/scATAC-seq/data/toy_data/peaks.bed" ;
std::string file_bam = "/local/groux/scATAC-seq/data/toy_data/sc_reads.bam" ;
std::string file_bai = "/local/groux/scATAC-seq/data/toy_data/sc_reads.bam.bai" ;

// GenomeRegion test suite
SUITE(GenomeRegion)
{
    // displays message
    TEST(message)
    {   std::cout << "Starting GenomicTools tests..." << std::endl ; }

    // tests vonstructor with value
    TEST(constructor_value)
    {   std::string chr = "chr1" ;
        int idx = 0 ;

        GenomeRegion r1(chr, idx, 0, 10) ;
        CHECK_EQUAL(chr, r1.chromosome) ;
        CHECK_EQUAL(0, r1.start) ;
        CHECK_EQUAL(10, r1.end) ;
        CHECK_EQUAL(10, r1.length) ;

        GenomeRegion r2(chr, idx, 1, 10) ;
        CHECK_EQUAL(chr, r2.chromosome) ;
        CHECK_EQUAL(1, r2.start) ;
        CHECK_EQUAL(10, r2.end) ;
        CHECK_EQUAL(9, r2.length) ;

        CHECK_THROW(GenomeRegion(chr, idx, -1,  10), std::invalid_argument) ;
        CHECK_THROW(GenomeRegion(chr, idx,  0, -10), std::invalid_argument) ;
    }


    // tests constructFragment factory function to create regions from bam
    /*
    TEST(test_contructFragment)
    {
        // expected content of bam file
        std::vector<GenomeRegion> regions ;
        regions.push_back(GenomeRegion("chr1",  400,   480)) ;
        regions.push_back(GenomeRegion("chr1",  470,   550)) ;
        regions.push_back(GenomeRegion("chr1",  560,   800)) ;
        regions.push_back(GenomeRegion("chr1",  560,   640)) ;
        regions.push_back(GenomeRegion("chr1",  610,   690)) ;
        regions.push_back(GenomeRegion("chr1",  670,   750)) ;
        regions.push_back(GenomeRegion("chr1",  730,   810)) ;
        regions.push_back(GenomeRegion("chr1",  770,   850)) ;
        regions.push_back(GenomeRegion("chr1",  950,   1150)) ;
        regions.push_back(GenomeRegion("chr1",  960,   1040)) ;
        regions.push_back(GenomeRegion("chr1",  1010,  1090)) ;
        regions.push_back(GenomeRegion("chr1",  1060,  1140)) ;
        regions.push_back(GenomeRegion("chr1",  1070,  1150)) ;
        regions.push_back(GenomeRegion("chr1",  1350,  1430)) ;
        regions.push_back(GenomeRegion("chr1",  1360,  1440)) ;
        regions.push_back(GenomeRegion("chr1",  1410,  1490)) ;
        regions.push_back(GenomeRegion("chr1",  1500,  1600)) ;
        regions.push_back(GenomeRegion("chr1",  1600,  1700)) ;

        regions.push_back(GenomeRegion("chr2",  400,   480)) ;
        regions.push_back(GenomeRegion("chr2",  470,   550)) ;
        regions.push_back(GenomeRegion("chr2",  560,   800)) ;
        regions.push_back(GenomeRegion("chr2",  560,   640)) ;
        regions.push_back(GenomeRegion("chr2",  610,   690)) ;
        regions.push_back(GenomeRegion("chr2",  670,   750)) ;
        regions.push_back(GenomeRegion("chr2",  730,   810)) ;
        regions.push_back(GenomeRegion("chr2",  770,   850)) ;
        regions.push_back(GenomeRegion("chr2",  950,   1150)) ;
        regions.push_back(GenomeRegion("chr2",  960,   1040)) ;
        regions.push_back(GenomeRegion("chr2",  1010,  1090)) ;
        regions.push_back(GenomeRegion("chr2",  1060,  1140)) ;
        regions.push_back(GenomeRegion("chr2",  1070,  1150)) ;
        regions.push_back(GenomeRegion("chr2",  1350,  1430)) ;
        regions.push_back(GenomeRegion("chr2",  1360,  1440)) ;
        regions.push_back(GenomeRegion("chr2",  1410,  1490)) ;
        regions.push_back(GenomeRegion("chr2",  1500,  1600)) ;
        regions.push_back(GenomeRegion("chr2",  1600,  1700)) ;

        seqan::BamAlignmentRecord record ;
        std::string bam_path = "src/Unittests/data/sc_reads.bam" ;

        // read file for fragments starting on + strand
        seqan::BamFileIn bam_file(bam_path.c_str()) ;
        // header
        seqan::BamHeader bam_header ;
        seqan::readHeader(bam_header, bam_file) ;
        for(size_t i=0; not seqan::atEnd(bam_file); i++)
        {   seqan::readRecord(record, bam_file) ;
            if(seqan::hasFlagFirst(record) and not seqan::hasFlagRC(record))
            {   std::cout << regions[i] << "   "
                          << GenomeRegion::constructFragment(record)
                          << std::endl ;
                CHECK_EQUAL(regions[i], GenomeRegion::constructFragment(record)) ; }
        }
        seqan::close(bam_file) ;

        // read file for fragments starting on - strand
        seqan::BamFileIn bam_file(bam_path.c_str()) ;
        // header
        seqan::BamHeader bam_header ;
        seqan::readHeader(bam_header, bam_file) ;
        for(size_t i=0; not seqan::atEnd(bam_file); i++)
        {   seqan::readRecord(record, bam_file) ;
            if(seqan::hasFlagFirst(record) and seqan::hasFlagRC(record))
            {   CHECK_EQUAL(regions[i], GenomeRegion::constructFragment(record)) ; }
        }
        seqan::close(bam_file) ;
    }
    */

    TEST(test_contructRead)
    {   // expected content of bam file
        std::list<GenomeRegion> regions_exp ;
        // chromosome 1 -> has index 0 in BAM file header
        regions_exp.push_back(GenomeRegion("chr1", 0, 400,  435)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 400,  435)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 445,  480)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 445,  480)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 470,  505)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 470,  505)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 515,  550)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 515,  550)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 560,  595)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 560,  595)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 560,  595)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 560,  595)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 605,  640)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 605,  640)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 610,  645)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 610,  645)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 655,  690)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 655,  690)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 670,  705)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 670,  705)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 715,  750)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 715,  750)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 730,  765)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 730,  765)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 765,  800)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 765,  800)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 770,  805)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 770,  805)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 775,  810)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 775,  810)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 815,  850)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 815,  850)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 950,  985)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 950,  985)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 960,  995)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 960,  995)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1005, 1040)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1005, 1040)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1010, 1045)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1010, 1045)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1055, 1090)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1055, 1090)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1060, 1095)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1060, 1095)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1070, 1105)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1070, 1105)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1105, 1140)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1105, 1140)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1115, 1150)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1115, 1150)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1115, 1150)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1115, 1150)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1350, 1385)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1350, 1385)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1360, 1395)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1360, 1395)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1395, 1430)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1395, 1430)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1405, 1440)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1405, 1440)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1410, 1445)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1410, 1445)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1455, 1490)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1455, 1490)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1500, 1535)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1500, 1535)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1565, 1600)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1565, 1600)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1600, 1635)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1600, 1635)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1665, 1700)) ;
        regions_exp.push_back(GenomeRegion("chr1", 0, 1665, 1700)) ;

        // chromosome 2 -> has index 1 in BAM file header
        regions_exp.push_back(GenomeRegion("chr2", 1, 400,  435)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 400,  435)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 445,  480)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 445,  480)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 470,  505)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 470,  505)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 515,  550)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 515,  550)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 560,  595)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 560,  595)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 560,  595)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 560,  595)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 605,  640)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 605,  640)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 610,  645)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 610,  645)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 655,  690)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 655,  690)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 670,  705)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 670,  705)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 715,  750)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 715,  750)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 730,  765)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 730,  765)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 765,  800)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 765,  800)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 770,  805)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 770,  805)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 775,  810)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 775,  810)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 815,  850)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 815,  850)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 950,  985)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 950,  985)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 960,  995)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 960,  995)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1005, 1040)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1005, 1040)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1010, 1045)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1010, 1045)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1055, 1090)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1055, 1090)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1060, 1095)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1060, 1095)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1070, 1105)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1070, 1105)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1105, 1140)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1105, 1140)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1115, 1150)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1115, 1150)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1115, 1150)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1115, 1150)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1350, 1385)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1350, 1385)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1360, 1395)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1360, 1395)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1395, 1430)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1395, 1430)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1405, 1440)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1405, 1440)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1410, 1445)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1410, 1445)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1455, 1490)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1455, 1490)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1500, 1535)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1500, 1535)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1565, 1600)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1565, 1600)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1600, 1635)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1600, 1635)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1665, 1700)) ;
        regions_exp.push_back(GenomeRegion("chr2", 1, 1665, 1700)) ;

        // open file
        seqan::BamFileIn bam_file ;
        if (!seqan::open(bam_file, file_bam.c_str()))
        {   char msg[4096] ;
            sprintf(msg, "ERROR: could not open input file %s", file_bam.c_str()) ;
        }

        // read file
        seqan::BamAlignmentRecord record ;
        seqan::BamHeader header ;
        seqan::readHeader(header, bam_file) ;
        std::list<GenomeRegion> regions_val ;
        while(not seqan::atEnd(bam_file))
        {   seqan::readRecord(record, bam_file) ;
            regions_val.push_back(GenomeRegion::constructRead(record, bam_file)) ;
        }
        seqan::close(bam_file) ;

        // compare
        CHECK_EQUAL(regions_exp.size(), regions_val.size()) ;
        auto iter_exp = regions_exp.begin() ;
        auto iter_val = regions_val.begin() ;
        while(iter_exp != regions_exp.end())
        {   CHECK_EQUAL(*iter_exp, *iter_val) ;
            iter_exp++ ;
            iter_val++ ;
        }
    }

    // tests the method to check overlaps
    TEST(overlap)
    {   GenomeRegion r1("chr1", 0, 20, 30) ; // reference
        GenomeRegion r2("chr1", 0, 20, 30) ; // same as reference
        GenomeRegion r3("chr1", 0,  0, 45) ; // totally contain reference
        GenomeRegion r4("chr1", 0,  0, 10) ; // no overlap, upstream reference
        GenomeRegion r5("chr1", 0, 15, 25) ; // partial overlap reference
        GenomeRegion r6("chr1", 0, 22, 29) ; // inside reference
        GenomeRegion r7("chr1", 0, 25, 35) ; // partial overlap reference
        GenomeRegion r8("chr1", 0, 35, 45) ; // no overlap, downstream reference
        GenomeRegion r9("chr2", 1, 20, 30) ; // diff chromosome

        // always check reciprocity
        CHECK_EQUAL(true,  r1 | r1) ;
        CHECK_EQUAL(true,  r1 | r2) ; CHECK_EQUAL(true,  r2 | r1) ;
        CHECK_EQUAL(true,  r1 | r3) ; CHECK_EQUAL(true,  r3 | r1) ;
        CHECK_EQUAL(false, r1 | r4) ; CHECK_EQUAL(false, r4 | r1) ;
        CHECK_EQUAL(true,  r1 | r5) ; CHECK_EQUAL(true,  r5 | r1) ;
        CHECK_EQUAL(true,  r1 | r6) ; CHECK_EQUAL(true,  r6 | r1) ;
        CHECK_EQUAL(true,  r1 | r7) ; CHECK_EQUAL(true,  r7 | r1) ;
        CHECK_EQUAL(false, r1 | r8) ; CHECK_EQUAL(false, r8 | r1) ;
        CHECK_EQUAL(false, r1 | r9) ; CHECK_EQUAL(false, r9 | r1) ;
    }

    // tests the methods to get overlap length
    TEST(overlap_len)
    {   GenomeRegion r1("chr1", 0, 10, 20) ; // reference
        GenomeRegion r2("chr1", 0, 10, 20) ; // same as reference
        GenomeRegion r3("chr1", 0,  0, 45) ; // totally contain reference
        GenomeRegion r4("chr2", 1, 10, 20) ; // diff chromosome

        // always check reciprocity
        CHECK_EQUAL(10, r1.overlap_len(r1)) ;
        CHECK_EQUAL(10, r1.overlap_len(r2)) ; CHECK_EQUAL(10, r1.overlap_len(r2)) ;
        CHECK_EQUAL(10, r1.overlap_len(r3)) ; CHECK_EQUAL(10, r1.overlap_len(r3)) ;
        CHECK_EQUAL(0,  r1.overlap_len(r4)) ; CHECK_EQUAL(0,  r1.overlap_len(r4)) ;

        // slide a smaller region along reference, from before to after
        std::vector<int> overlaps = {0,0,1,2,3,4,4,4,4,4,4,4,3,2,1,0,0,0} ;
        int len = 4 ;
        for(int i=0, start=5; start<23; i++, start++)
        {   int end = start + len ;
            GenomeRegion s1("chr1", 0, start, end) ;
            CHECK_EQUAL(overlaps[i], r1.overlap_len(s1)) ;
            CHECK_EQUAL(overlaps[i], s1.overlap_len(r1)) ;
        }
    }

    // tests the is upstream and is downstream operators
    TEST(upstream_downstream)
    {   GenomeRegion r1("chr1", 0, 10, 20) ; // reference
        GenomeRegion r2("chr1", 0, 10, 20) ; // same as reference
        GenomeRegion r3("chr1", 0,  0, 45) ; // totally contain reference
        GenomeRegion r4("chr2", 1, 10, 20) ; // diff chromosome (downstream has 0 < 1)

        // always check reciprocity
        CHECK_EQUAL(false, r1 < r1) ; CHECK_EQUAL(false, r1 > r1) ;
        CHECK_EQUAL(false, r1 < r2) ; CHECK_EQUAL(false, r1 > r2) ;

        CHECK_EQUAL(false, r1 < r3) ; CHECK_EQUAL(false, r1 < r3) ;
        CHECK_EQUAL(false, r3 < r1) ; CHECK_EQUAL(false, r3 < r1) ;

        // not on the same chromosome -> depends on the index value
        CHECK_EQUAL(r1 < r4, true)  ; CHECK_EQUAL(r1 > r4, false) ;
        CHECK_EQUAL(r4 < r1, false) ; CHECK_EQUAL(r4 > r1, true) ;

        // slide a smaller region along reference, from before to after
        std::vector<bool> s1_upstream   = {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} ; // s1 < r1
        std::vector<bool> r1_downstream = {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} ; // r1 > s1
        std::vector<bool> s1_downstream = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1} ; // s1 > r1
        std::vector<bool> r1_upstream   = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1} ; // r1 < s1
        int len = 4 ;
        for(int i=0, start=5; start<23; i++, start++)
        {   // the sliding one
            int end = start + len ;
            GenomeRegion s1("chr1", 0, start, end) ;

            CHECK_EQUAL(s1_upstream[i],   s1 < r1) ;
            CHECK_EQUAL(r1_downstream[i], r1 > s1) ;

            CHECK_EQUAL(r1_upstream[i],   r1 < s1) ;
            CHECK_EQUAL(s1_downstream[i], s1 > r1) ;
        }
    }
}


// CorrelationMatrixCreator test suite
SUITE(CorrelationMatrixCreator)
{
    // displays message
    TEST(message)
    {   std::cout << "Starting CorrelationMatrixCreator tests..." << std::endl ; }

    // tests matrix creation with full fragments
    TEST(create_matrix_fragment)
    {   CorrelationMatrixCreator creator(file_bed,
                                            file_bam,
                                            file_bai,
                                            -500,
                                             500,
                                             100,
                                             CorrelationMatrixCreator::FRAGMENT) ;
        Matrix2D<int> m_val = creator.create_matrix() ;
        Matrix2D<int> m_exp(2, 9, 0) ;
        m_exp(0,0) = 420 ; m_exp(0,1) = 480 ; m_exp(0,2) = 380 ;
        m_exp(0,3) =   0 ; m_exp(0,4) = 440 ; m_exp(0,5) = 600 ;
        m_exp(0,6) =   0 ; m_exp(0,7) =   0 ; m_exp(0,8) = 400 ;

        m_exp(1,0) = 420 ; m_exp(1,1) = 480 ; m_exp(1,2) = 380 ;
        m_exp(1,3) =   0 ; m_exp(1,4) = 440 ; m_exp(1,5) = 600 ;
        m_exp(1,6) =   0 ; m_exp(1,7) =   0 ; m_exp(1,8) = 400 ;

        CHECK_EQUAL(m_exp.get_nrow(), m_val.get_nrow()) ;
        CHECK_EQUAL(m_exp.get_ncol(), m_val.get_ncol()) ;

        for(size_t i=0; i<m_exp.get_nrow(); i++)
        {   for(size_t j=0; j<m_exp.get_ncol(); j++)
            {   CHECK_EQUAL(m_exp(i,j), m_val(i,j)) ; }
        }
    }

    // tests matrix creation with fragment centers
    TEST(create_matrix_fragment_center)
    {   CorrelationMatrixCreator creator(file_bed,
                                            file_bam,
                                            file_bai,
                                            -500,
                                             500,
                                             100,
                                            CorrelationMatrixCreator::FRAGMENT_CENTER) ;
        Matrix2D<int> m_val = creator.create_matrix() ;
        Matrix2D<int> m_exp(2, 9, 0) ;
        m_exp(0,0) = 2 ; m_exp(0,1) = 6 ; m_exp(0,2) = 4 ;
        m_exp(0,3) = 0 ; m_exp(0,4) = 2 ; m_exp(0,5) = 8 ;
        m_exp(0,6) = 0 ; m_exp(0,7) = 0 ; m_exp(0,8) = 4 ;

        m_exp(1,0) = 2 ; m_exp(1,1) = 6 ; m_exp(1,2) = 4 ;
        m_exp(1,3) = 0 ; m_exp(1,4) = 2 ; m_exp(1,5) = 8 ;
        m_exp(1,6) = 0 ; m_exp(1,7) = 0 ; m_exp(1,8) = 4 ;

        CHECK_EQUAL(m_exp.get_nrow(), m_val.get_nrow()) ;
        CHECK_EQUAL(m_exp.get_ncol(), m_val.get_ncol()) ;

        for(size_t i=0; i<m_exp.get_nrow(); i++)
        {   for(size_t j=0; j<m_exp.get_ncol(); j++)
            {   CHECK_EQUAL(m_exp(i,j), m_val(i,j)) ; }
        }
    }

    // tests matrix creation with reads
    TEST(create_matrix_read)
    {   CorrelationMatrixCreator creator(file_bed,
                                            file_bam,
                                            file_bai,
                                            -500,
                                             500,
                                             100,
                                            CorrelationMatrixCreator::READ) ;
        Matrix2D<int> m_val = creator.create_matrix() ;
        Matrix2D<int> m_exp(2, 9, 0) ;
        m_exp(0,0) = 280 ; m_exp(0,1) = 250 ; m_exp(0,2) = 310 ;
        m_exp(0,3) =   0 ; m_exp(0,4) = 280 ; m_exp(0,5) = 420 ;
        m_exp(0,6) =   0 ; m_exp(0,7) =   0 ; m_exp(0,8) = 350 ;

        m_exp(1,0) = 280 ; m_exp(1,1) = 250 ; m_exp(1,2) = 310 ;
        m_exp(1,3) =   0 ; m_exp(1,4) = 280 ; m_exp(1,5) = 420 ;
        m_exp(1,6) =   0 ; m_exp(1,7) =   0 ; m_exp(1,8) = 350 ;

        CHECK_EQUAL(m_exp.get_nrow(), m_val.get_nrow()) ;
        CHECK_EQUAL(m_exp.get_ncol(), m_val.get_ncol()) ;

        for(size_t i=0; i<m_exp.get_nrow(); i++)
        {   for(size_t j=0; j<m_exp.get_ncol(); j++)
            {   CHECK_EQUAL(m_exp(i,j), m_val(i,j)) ; }
        }
    }


    // tests matrix creation with ATAC-seq reads
    TEST(create_matrix_read_atac)
    {   CorrelationMatrixCreator creator(file_bed,
                                            file_bam,
                                            file_bai,
                                            -500,
                                             500,
                                             100,
                                            CorrelationMatrixCreator::READ_ATAC) ;
        Matrix2D<int> m_val = creator.create_matrix() ;
        Matrix2D<int> m_exp(2, 9, 0) ;
        m_exp(0,0) = 8 ; m_exp(0,1) = 8 ; m_exp(0,2) =  8 ;
        m_exp(0,3) = 0 ; m_exp(0,4) = 8 ; m_exp(0,5) = 12 ;
        m_exp(0,6) = 0 ; m_exp(0,7) = 0 ; m_exp(0,8) = 10 ;

        m_exp(1,0) = 8 ; m_exp(1,1) = 8 ; m_exp(1,2) =  8 ;
        m_exp(1,3) = 0 ; m_exp(1,4) = 8 ; m_exp(1,5) = 12 ;
        m_exp(1,6) = 0 ; m_exp(1,7) = 0 ; m_exp(1,8) = 10 ;

        CHECK_EQUAL(m_exp.get_nrow(), m_val.get_nrow()) ;
        CHECK_EQUAL(m_exp.get_ncol(), m_val.get_ncol()) ;

        for(size_t i=0; i<m_exp.get_nrow(); i++)
        {   for(size_t j=0; j<m_exp.get_ncol(); j++)
            {   CHECK_EQUAL(m_exp(i,j), m_val(i,j)) ; }
        }
    }
}
