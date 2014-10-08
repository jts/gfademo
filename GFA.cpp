//-------------------------------------------------------------------------------
// 
// GFA - Graphical Fragment Assembly format parser 
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
#include <fstream>
#include <iostream>
#include <assert.h>

//
//
//
void GFA::readFile(const std::string& filename,
                   SegmentVector& segments,
                   LinkVector& links,
                   ContainmentVector& containments)
{
    std::ifstream input(filename.c_str());
    
    char record_type;
    while(input >> record_type) {
        switch(record_type)
        {
            case 'H':
                parseHeader(input);
                break;
            case 'S':
                segments.push_back(parseSegment(input));
                break;
            case 'L':
                links.push_back(parseLink(input));
                break;
            case 'C':
                containments.push_back(parseContainment(input));
                break;
            default:
                printf("Unknown tag: %c\n", record_type);
                assert(false);
                break;
        }
    }
}

GFA::Header GFA::parseHeader(std::istream& input)
{
    GFA::Header header;
    getline(input, header.text);
    return header;
}

GFA::Segment GFA::parseSegment(std::istream& input)
{
    GFA::Segment segment;
    input >> segment.id;
    input >> segment.sequence;

    // Parse any tags but discard them for now
    std::string tags;
    getline(input, tags);

    return segment;
}

GFA::Link GFA::parseLink(std::istream& input)
{
    GFA::Link link;
    input >> link.id[0];
    input >> link.orientation[0];
    input >> link.id[1];
    input >> link.orientation[1];
    input >> link.cigar;
    return link;
}

GFA::Containment GFA::parseContainment(std::istream& input)
{
    GFA::Containment containment;
    input >> containment.id[0];
    input >> containment.orientation[0];
    input >> containment.id[1];
    input >> containment.orientation[2];
    input >> containment.offset;
    input >> containment.cigar;
    return containment;
}
