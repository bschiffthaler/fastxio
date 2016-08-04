// ============================================================================
// b2zstream, C++ iostream classes wrapping the bzlib compression library.
// Copyright (C) 2001 Deepak Bandyopadhyay, Lutz Kettner, Bastian Schiffthaler
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : bz2stream.h
// Revision      : $Revision: 1
// Revision_date : $Date: 2016/14/07
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner, Bastian Schiffthaler
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library", copied almost verbatim from Deepak Bandyopadhyay's
// and Lutz Kettner's gzstream library.
// ============================================================================

#ifndef BZ2STREAM_H
#define BZ2STREAM_H 1

// standard C++ with new header file names and std:: namespace
#include <iostream>
#include <fstream>
#include <bzlib.h>

#ifdef BZ2STREAM_NAMESPACE
namespace BZ2STREAM_NAMESPACE {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement bz2stream. See below for user classes.
// ----------------------------------------------------------------------------

class bz2streambuf : public std::streambuf {
private:
    static const int bufferSize = 47+256;    // size of data buff
    // totals 512 bytes under g++ for igzstream at the end.

    BZFILE*          file;               // file handle for compressed file
    char             buffer[bufferSize]; // data buffer
    char             opened;             // open/close state of stream
    int              mode;               // I/O mode

    int flush_buffer();
public:
    bz2streambuf() : opened(0) {
        setp( buffer, buffer + (bufferSize-1));
        setg( buffer + 4,     // beginning of putback area
              buffer + 4,     // read position
              buffer + 4);    // end position      
        // ASSERT: both input & output capabilities will not be used together
    }
    int is_open() { return opened; }
    bz2streambuf* open( const char* name, int open_mode);
    bz2streambuf* close();
    ~bz2streambuf() { close(); }
    
    virtual int     overflow( int c = EOF);
    virtual int     underflow();
    virtual int     sync();
};

class bz2streambase : virtual public std::ios {
protected:
    bz2streambuf buf;
public:
    bz2streambase() { init(&buf); }
    bz2streambase( const char* name, int open_mode);
    ~bz2streambase();
    void open( const char* name, int open_mode);
    void close();
    bz2streambuf* rdbuf() { return &buf; }
};

// ----------------------------------------------------------------------------
// User classes. Use ibz2stream and obz2stream analogously to ifstream and
// ofstream respectively. They read and write files based on the BZ2_bz* 
// function interface of the bzlib. Files are compatible with bzip2 compression.
// ----------------------------------------------------------------------------

class ibz2stream : public bz2streambase, public std::istream {
public:
    ibz2stream() : std::istream( &buf) {} 
    ibz2stream( const char* name, int open_mode = std::ios::in)
        : bz2streambase( name, open_mode), std::istream( &buf) {}  
    bz2streambuf* rdbuf() { return bz2streambase::rdbuf(); }
    void open( const char* name, int open_mode = std::ios::in) {
        bz2streambase::open( name, open_mode);
    }
};

class obz2stream : public bz2streambase, public std::ostream {
public:
    obz2stream() : std::ostream( &buf) {}
    obz2stream( const char* name, int mode = std::ios::out)
        : bz2streambase( name, mode), std::ostream( &buf) {}  
    bz2streambuf* rdbuf() { return bz2streambase::rdbuf(); }
    void open( const char* name, int open_mode = std::ios::out) {
        bz2streambase::open( name, open_mode);
    }
};

#ifdef BZ2STREAM_NAMESPACE
} // namespace GZSTREAM_NAMESPACE
#endif

#endif // BZ2STREAM_H
// ============================================================================
// EOF //

