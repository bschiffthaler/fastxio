// ============================================================================
// bz2stream, C++ iostream classes wrapping the bzlib compression library.
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner, Bastian Schiffthaler
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
// File          : bz2stream.cpp
// Revision      : $Revision: 1 $
// Revision_date : $Date: 2016/14/07
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner, Bastian Schiffthaler
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library", copied almost verbatim from Deepak Bandyopadhyay's
// and Lutz Kettner's gzstream library.
// ============================================================================

#include <bzlib.h>
#include <bz2stream.h>
#include <iostream>
#include <string.h>  // for memcpy

#ifdef BZ2STREAM_NAMESPACE
namespace BZ2STREAM_NAMESPACE {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement bz2stream. See header file for user classes.
// ----------------------------------------------------------------------------

// --------------------------------------
// class bz2streambuf:
// --------------------------------------

bz2streambuf* bz2streambuf::open( const char* name, int open_mode) {
    if ( is_open())
        return (bz2streambuf*)0;
    mode = open_mode;
    // no append nor read/write mode
    if ((mode & std::ios::ate) || (mode & std::ios::app)
        || ((mode & std::ios::in) && (mode & std::ios::out)))
        return (bz2streambuf*)0;
    char  fmode[10];
    char* fmodeptr = fmode;
    if ( mode & std::ios::in)
        *fmodeptr++ = 'r';
    else if ( mode & std::ios::out)
        *fmodeptr++ = 'w';
    *fmodeptr++ = 'b';
    *fmodeptr = '\0';
    file = BZ2_bzopen( name, fmode);
    if (file == 0)
        return (bz2streambuf*)0;
    opened = 1;
    return this;
}

bz2streambuf * bz2streambuf::close() {
    if ( is_open()) {
        sync();
        opened = 0;
        if ( BZ2_bzflush( file) == BZ_OK) {
	    BZ2_bzclose( file);
            return this;
	}
    }
    return (bz2streambuf*)0;
}

int bz2streambuf::underflow() { // used for input buffer only
    if ( gptr() && ( gptr() < egptr()))
        return * reinterpret_cast<unsigned char *>( gptr());

    if ( ! (mode & std::ios::in) || ! opened)
        return EOF;
    // Josuttis' implementation of inbuf
    int n_putback = gptr() - eback();
    if ( n_putback > 4)
        n_putback = 4;
    memcpy( buffer + (4 - n_putback), gptr() - n_putback, n_putback);

    int num = BZ2_bzread( file, buffer+4, bufferSize-4);
    if (num <= 0) // ERROR or EOF
        return EOF;

    // reset buffer pointers
    setg( buffer + (4 - n_putback),   // beginning of putback area
          buffer + 4,                 // read position
          buffer + 4 + num);          // end of buffer

    // return next character
    return * reinterpret_cast<unsigned char *>( gptr());    
}

int bz2streambuf::flush_buffer() {
    // Separate the writing of the buffer from overflow() and
    // sync() operation.
    int w = pptr() - pbase();
    if ( BZ2_bzwrite( file, pbase(), w) != w)
        return EOF;
    pbump( -w);
    return w;
}

int bz2streambuf::overflow( int c) { // used for output buffer only
    if ( ! ( mode & std::ios::out) || ! opened)
        return EOF;
    if (c != EOF) {
        *pptr() = c;
        pbump(1);
    }
    if ( flush_buffer() == EOF)
        return EOF;
    return c;
}

int bz2streambuf::sync() {
    // Changed to use flush_buffer() instead of overflow( EOF)
    // which caused improper behavior with std::endl and flush(),
    // bug reported by Vincent Ricard.
    if ( pptr() && pptr() > pbase()) {
        if ( flush_buffer() == EOF)
            return -1;
    }
    return 0;
}

// --------------------------------------
// class bz2streambase:
// --------------------------------------

bz2streambase::bz2streambase( const char* name, int mode) {
    init( &buf);
    open( name, mode);
}

bz2streambase::~bz2streambase() {
    buf.close();
}

void bz2streambase::open( const char* name, int open_mode) {
    if ( ! buf.open( name, open_mode))
        clear( rdstate() | std::ios::badbit);
}

void bz2streambase::close() {
    if ( buf.is_open())
        if ( ! buf.close())
            clear( rdstate() | std::ios::badbit);
}

#ifdef BZ2STREAM_NAMESPACE
} // namespace BZ2STREAM_NAMESPACE
#endif

// ============================================================================
// EOF //
