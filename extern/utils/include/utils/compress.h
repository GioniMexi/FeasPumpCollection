/**
 * @file compress.h
 * @brief GZIP compression utilities
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2013
 */

#ifndef COMPRESS_H
#define COMPRESS_H

#include <iostream>
#include <string>
#include <cstdio>
#include <zlib.h>

namespace dominiqs {

/**
 * Compress the given file @param inPath, storing the result in file @param outPath and using mode @param mode.
 * Remove the original if @param removeOriginal is true.
 */
void gzipCompress(const std::string& inPath, const std::string& outPath,
	const std::string& mode = "wb6", bool removeOriginal = true);

/**
 * Uncompress the given file @param inPath, storing the result in file @param outPath.
 * Remove the original if @param removeOriginal is true.
 */
void gzipUncompress(const std::string& inPath, const std::string& outPath, bool removeOriginal = true);


/**
 * Compressed Stream Buffer
 */

class gzstreambuf : public std::streambuf
{
private:
	static const int GZ_BUFFER_SIZE = 47+256;    // size of data buff
	// totals 512 bytes under g++ for igzstream at the end.
	gzFile file; // file handle for compressed file
	char buffer[GZ_BUFFER_SIZE]; // data buffer
	bool opened; // open/close state of stream
	int mode;  // GZIP I/O mode
	int flush_buffer();
public:
	gzstreambuf() : opened(false)
	{
		setp(buffer, buffer + (GZ_BUFFER_SIZE - 1));
		setg(buffer + 4, // beginning of putback area
			buffer + 4, // read position
			buffer + 4); // end position
		// ASSERT: both input & output capabilities will not be used together
	}
	int is_open() { return opened; }
	gzstreambuf* open(const char* name, int open_mode);
	gzstreambuf* close();
	~gzstreambuf() { close(); }

	virtual int overflow(int c = EOF);
	virtual int underflow();
	virtual int sync();
};

/**
 * Base Class for GZStream classes
 */
class gzstreambase : virtual public std::ios
{
protected:
	gzstreambuf buf;
public:
	gzstreambase() { init(&buf); }
	gzstreambase(const char* name, int mode)  { init(&buf); open(name, mode); }
	~gzstreambase() { buf.close(); }
	void open(const char* name, int mode)
	{
		if (!buf.open(name, mode)) clear(rdstate() | std::ios::badbit);
	}
	void close()
	{
		if (buf.is_open() && !buf.close()) clear(rdstate() | std::ios::badbit);
	}
	gzstreambuf* rdbuf() { return &buf; }
};

/**
 * Compressed Input File Stream
 */
class igzstream : public gzstreambase, public std::istream
{
public:
	igzstream() : std::istream(&buf) {} 
	igzstream(const char* name, int open_mode = std::ios::in) : gzstreambase(name, open_mode), std::istream(&buf) {}  
	gzstreambuf* rdbuf() { return gzstreambase::rdbuf(); }
	void open(const char* name, int open_mode = std::ios::in) { gzstreambase::open(name, open_mode); }
};

/**
 * Compressed Output File Stream
 */
class ogzstream : public gzstreambase, public std::ostream
{
public:
	ogzstream() : std::ostream(&buf) {}
	ogzstream(const char* name, int mode = std::ios::out) : gzstreambase(name, mode), std::ostream(&buf) {}  
	gzstreambuf* rdbuf() { return gzstreambase::rdbuf(); }
	void open(const char* name, int open_mode = std::ios::out) { gzstreambase::open(name, open_mode); }
};

} // namespace dominiqs

#endif /* end of include guard: COMPRESS_H */
