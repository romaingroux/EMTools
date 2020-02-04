#include <GenomeRegion.hpp>

#include <string>
#include <cmath>          // abs()
#include <stdexcept>      // std::invalid_argument
#include <seqan/bam_io.h>


GenomeRegion GenomeRegion::constructRead(const seqan::BamAlignmentRecord& record,
                                         const seqan::BamFileIn& bam_file)
{
    GenomeRegion read ;

    read.chromosome =
                seqan::toCString(
                seqan::getContigName(record, bam_file)) ;
    read.chromosome_idx = record.rID ;

    read.start  = record.beginPos ;
    read.length = seqan::endPosition(record.seq) ;
    read.end    = read.start + read.length ;

    if(read.start < 0 or read.end < 0)
    {   char msg[4096] ;
        sprintf(msg, "Error! invalide coordinate (<0) : [%s/%d %d %d)]",
                read.chromosome.c_str(), read.chromosome_idx,
                read.start, read.end) ;
        throw std::invalid_argument(msg) ;
    }
    else if(read.start >= read.end)
    {   char msg[4096] ;
        sprintf(msg, "Error! start >= end : [%s/%d %d %d)]",
                read.chromosome.c_str(), read.chromosome_idx,
                read.start, read.end) ;
        throw std::invalid_argument(msg) ;
    }
    return read ;
}

GenomeRegion GenomeRegion::constructReadATAC(const seqan::BamAlignmentRecord& record,
                                             const seqan::BamFileIn& bam_file)
{   GenomeRegion read = GenomeRegion::constructRead(record, bam_file);
    if(not seqan::hasFlagRC(record))
    {   read.start  += 4 ;
        read.end    = read.start + 1 ;
        read.length = 1 ;
    }
    else
    {   read.start  = read.end - 1 - 5  ;
        read.end    = read.start + 1 ;
        read.length = 1 ;
    }
    return read ;
}

GenomeRegion GenomeRegion::constructReadEdge(const seqan::BamAlignmentRecord& record,
                                             const seqan::BamFileIn& bam_file)
{
    GenomeRegion read = GenomeRegion::constructRead(record, bam_file);
    if(not seqan::hasFlagRC(record))
    {   read.end    = read.start + 1 ;
        read.length = 1 ;
    }
    else
    {   read.start  = read.end - 1  ;
        read.length = 1 ;
    }
    return read ;
}

GenomeRegion GenomeRegion::constructFragment(const seqan::BamAlignmentRecord& record,
                                             const seqan::BamFileIn& bam_file)
{   GenomeRegion frag ;

    frag.chromosome =
                seqan::toCString(
                seqan::getContigName(record, bam_file)) ;
    frag.chromosome_idx = record.rID ;


    // read is on + strand
    //  record
    // |----->    <-----|
    if(not seqan::hasFlagRC(record))
    {   frag.start  = record.beginPos ;
        frag.length = record.tLen ;
        frag.end    = frag.start + frag.length ;
    }
    // read is on - strand
    //            record
    // |----->    <-----|
    else
    {   // frag.end    = seqan::endPosition(record.seq) + 1 ;
        // frag.length = abs(record.tLen) ;
        // frag.start  = frag.end - frag.length ;
        frag.length = abs(record.tLen) ;
        frag.start  = record.pNext ;
        frag.end    = frag.start + frag.length ;
    }

    if(frag.start < 0 or frag.end < 0)
    {   char msg[4096] ;
        sprintf(msg, "Error! invalide coordinate (<0) : [%s/%d %d %d)]",
                frag.chromosome.c_str(), frag.chromosome_idx,
                frag.start, frag.end) ;
        throw std::invalid_argument(msg) ;
    }
    else if(frag.start >= frag.end)
    {   char msg[4096] ;
        sprintf(msg, "Error! start >= end : [%s/%d %d %d)]",
                frag.chromosome.c_str(), frag.chromosome_idx,
                frag.start, frag.end) ;
        throw std::invalid_argument(msg) ;
    }
    return frag ;
}

GenomeRegion GenomeRegion::constructFragmentCenter(const seqan::BamAlignmentRecord& record,
                                                   const seqan::BamFileIn& bam_file)
{   GenomeRegion frag = GenomeRegion::constructFragment(record, bam_file) ;
    int mid = frag.start + (frag.length / 2) ;
    frag.start  = mid ;
    frag.end    = mid + 1;
    frag.length = 1 ;
    return frag ;
}

GenomeRegion::GenomeRegion(const GenomeRegion& other)
    : chromosome(other.chromosome),
      chromosome_idx(other.chromosome_idx),
      start(other.start),
      end(other.end),
      length(other.length)
{}

GenomeRegion::GenomeRegion(const std::string& chromosome,
                           int chromosome_idx,
                           int start,
                           int end)
    : chromosome(chromosome),
      chromosome_idx(chromosome_idx),
      start(start),
      end(end),
      length(end - start)
{   if(this->start < 0 or this->end < 0)
    {   char msg[4096] ;
        sprintf(msg, "Error! invalide coordinate (<0) : [%s/%d %d %d)]",
                this->chromosome.c_str(), this->chromosome_idx,
                this->start, this->end) ;
        throw std::invalid_argument(msg) ;
    }
    else if(start >= end)
    {   char msg[4096] ;
        sprintf(msg, "Error! start >= end : [%s/%d %d %d)]",
                this->chromosome.c_str(), this->chromosome_idx,
                this->start, this->end) ;
        throw std::invalid_argument(msg) ;
    }
}

int GenomeRegion::overlap_len(const GenomeRegion& other) const
{   int len = 0 ;
    if((*this) | other)
    {   // this is contained in other or overlap perfectly other
        if(this->start >= other.start and this->end <= other.end)
        {   len = this->length ; }
        // start of this overlaps end other
        else if((other.start < this->start) and (other.end-1 >= this->start))
        {   len = other.end - this->start ; }
        // other contained in this (perect overlap is handled in first case)
        else if(other.start >= this->start and other.end <= this->end)
        {   len = other.length ; }
        // end of this overlaps start of other (only case left)
        else
        {   len = this->end - other.start ; }
    }
    return len ;
}

GenomeRegion& GenomeRegion::operator = (const GenomeRegion& rhs)
{   if(this == &rhs)
    {   return *this ; }
    this->start          = rhs.start ;
    this->end            = rhs.end ;
    this->length         = rhs.length ;
    this->chromosome     = rhs.chromosome ;
    this->chromosome_idx = rhs.chromosome_idx ;
    return *this ;
}

GenomeRegion& GenomeRegion::operator = (GenomeRegion&& rhs)
{   if(this == &rhs)
    {   return *this ; }
    this->start          = rhs.start ;
    this->end            = rhs.end ;
    this->length         = rhs.length ;
    this->chromosome     = rhs.chromosome ;
    this->chromosome_idx = rhs.chromosome_idx ;
    return *this ;
}

bool GenomeRegion::operator == (const GenomeRegion& rhs) const
{   if(this == &rhs)
    {   return true ; }
    if(this->chromosome_idx == rhs.chromosome_idx and
       this->start          == rhs.start and
       this->end            == rhs.end and
       this->length         == rhs.length)
    {   return true ; }
    return false ;
}

bool GenomeRegion::operator | (const GenomeRegion& rhs) const
{

    if((this->chromosome_idx != rhs.chromosome_idx) or // on diff chromosomes
       (rhs.end-1 < this->start) or                    // rhs upstream this
       (this->end-1 < rhs.start))                      // rhs downstream this
    {   return false ; }
    return true ;
}

bool GenomeRegion::operator < (const GenomeRegion& rhs) const
{
    if(this->chromosome_idx < rhs.chromosome_idx)
    {   return true ; }
    else if((this->chromosome_idx == rhs.chromosome_idx) and
            (this->end-1 < rhs.start))
    { return true ; }
    return false ;
}

bool GenomeRegion::operator > (const GenomeRegion& rhs) const
{   if(this->chromosome_idx > rhs.chromosome_idx)
    {   return true ; }
    if((this->chromosome_idx == rhs.chromosome_idx) and
       (rhs.end-1 < this->start))
    {   return true ; }
    return false ;
}

std::ostream& operator << (std::ostream& stream, const GenomeRegion& region)
{   stream << "( "
           << region.chromosome << "/"
           << region.chromosome_idx << " "
           << region.start << " "
           << region.end << " "
           << region.length << " ) " ;
    return stream ;
}
