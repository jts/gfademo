//-------------------------------------------------------------------------------
// 
// gfademo - small demonstration of the graphical fragment assembly format
//
// Copyright (C) 2014 Jared Simpson (jared.simpson@gmail.com)
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// ------------------------------------------------------------------------------
#include "GFA.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <map>

// Complement a base
inline char complement(char base)
{
    switch(base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        default:
            assert(false && "Unknown base!");
            return 'N';
    }
}

// Reverse complement a sequence
std::string reverseComplement(const std::string& seq)
{
    std::string out(seq.length(), 'A');
    size_t last_pos = seq.length() - 1;
    for(int i = last_pos; i >= 0; --i)
        out[last_pos - i] = complement(seq[i]);
    return out;
}

// Calculate how many bases of each sequence are consumed by the
// cigar string
void calculateCigarConsumption(const std::string& cigar,
                               int* consumed_1,
                               int* consumed_2)
{
    *consumed_1 = 0;
    *consumed_2 = 0;
    
    // parse cigar
    int len;
    char op;
    std::stringstream parser(cigar);
    while(parser >> len >> op) {
        assert(std::string("MIDS").find(op) != std::string::npos);
        if(op == 'M' || op == 'D' || op == 'S')
            *consumed_1 += len;
        if(op == 'M' || op == 'I' || op == 'S')
            *consumed_2 += len;
    }
}

void printAlignment(const std::string& s0, const std::string& s1,
                    int offset, const std::string& cigar)
{
    int curr_0 = offset;
    int curr_1 = 0;
    assert(curr_0 >= 0);

    // Parse CIGAR into an alignment
    std::string m0;
    std::string m1;

    // Add leading unmatched sequence
    m0.append(s0.substr(0, curr_0));
    m1.append(curr_0, '.');

    int length;
    char op;
    std::stringstream parser(cigar);

    while(parser >> length >> op) {

        switch(op) {
            case 'D':
                m0.append(s0.substr(curr_0, length));
                m1.append(length, '-');
                curr_0 += length;
                break;
            case 'I':
                m0.append(length, '-');
                m1.append(s1.substr(curr_1, length));
                curr_1 += length;
                break;
            case 'M':
                m0.append(s0.substr(curr_0, length));
                m1.append(s1.substr(curr_1, length));
                curr_0 += length;
                curr_1 += length;
                break;
        } 
    }

    // Add trailing unmatched sequence
    m0.append(s1.length() - curr_1, '.');
    m1.append(s1.substr(curr_1));

    int matches = 0;
    int mismatches = 0;
    int gaps = 0;
    for(size_t i = 0; i < m0.size(); ++i) {
        if(m0[i] == '.' || m1[i] == '.')
            continue;

        bool is_gap = (m0[i] == '-' || m1[i] == '-');
        gaps += is_gap;
        matches += (m0[i] == m1[i]);
        mismatches += (!is_gap && m0[i] != m1[i]);
    }

    double pi = 0.0f;
    if(matches > 0)
        pi = (double)matches / (gaps + matches + mismatches);

    // Print M0 and M1 in fixed-width lines
    size_t ms = m0.size();
    size_t fw = 120;

    for(size_t i = 0; i < ms; i += fw) {
        size_t end = std::min(ms, i + fw);
        std::string sub0 = m0.substr(i, end - i);
        std::string sub1 = m1.substr(i, end - i);
        printf("M0: %s\n", sub0.c_str());
        printf("M1: %s\n", sub1.c_str());
    }
    
    printf("Matches: %d Mismatches: %d Gaps: %d Identity: %.3lf\n", matches, mismatches, gaps, pi);
}

int main(int argc, char** argv)
{
    if(argc != 2) {
        fprintf(stderr, "Usage: ./gfademo file.gfa\n");
        exit(EXIT_FAILURE);
    }

    GFA::SegmentVector segments;
    GFA::LinkVector links;
    GFA::ContainmentVector containments;
    GFA::readFile(argv[1],
                  segments,
                  links,
                  containments);

    // Build a map from segment id -> *segment
    std::map<std::string, GFA::Segment*> segment_map;
    for(size_t i = 0; i < segments.size(); ++i)
        segment_map.insert(std::make_pair(segments[i].id, &segments[i]));

    // Print links
    for(size_t li = 0; li < links.size(); ++li) {

        GFA::Link link = links[li];
        printf("====== Link %s %c %s %c %s =======\n", 
               link.id[0].c_str(), link.orientation[0], 
               link.id[1].c_str(), link.orientation[1],
               link.cigar.c_str());

        // Determine how long the alignments are
        int n_aligned_0 = 0;
        int n_aligned_1 = 0;
        
        calculateCigarConsumption(link.cigar, &n_aligned_0, &n_aligned_1);

        GFA::Segment* segment_0 = segment_map[link.id[0]];
        GFA::Segment* segment_1 = segment_map[link.id[1]];
        assert(segment_0 && segment_1);

        std::string s0 = segment_0->sequence;
        std::string s1 = segment_1->sequence;
        
        assert(n_aligned_0 < s0.length());
        assert(n_aligned_1 < s1.length());

        if(link.orientation[0] == '-')
            s0 = reverseComplement(s0);
        
        if(link.orientation[1] == '-')
            s1 = reverseComplement(s1);
        
        printAlignment(s0, s1, s0.length() - n_aligned_0, link.cigar);
        printf("\n");
    }

    // Print Containments
    // TODO: refactor links and containments into common datastructure
    for(size_t ci = 0; ci < containments.size(); ++ci) {

        GFA::Containment containment = containments[ci];
        printf("====== Containment %s %c %s %c %d %s =======\n", 
               containment.id[0].c_str(), containment.orientation[0], 
               containment.id[1].c_str(), containment.orientation[1],
               containment.offset,
               containment.cigar.c_str());

        GFA::Segment* segment_0 = segment_map[containment.id[0]];
        GFA::Segment* segment_1 = segment_map[containment.id[1]];
        assert(segment_0 && segment_1);

        std::string s0 = segment_0->sequence;
        std::string s1 = segment_1->sequence;

        if(containment.orientation[0] == '-')
            s0 = reverseComplement(s0);
        
        if(containment.orientation[1] == '-')
            s1 = reverseComplement(s1);
        
        printAlignment(s0, s1, containment.offset, containment.cigar);
        printf("\n");
    }    
}
