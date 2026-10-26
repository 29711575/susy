#ifndef KITAEV_SSH_CHECKPOINT_IO_H
#define KITAEV_SSH_CHECKPOINT_IO_H

#include <cstdint>
#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <locale>

namespace kitaev_ssh {
inline void checkpointRequire(bool ok,const char *message) {
    if(!ok)throw std::invalid_argument(std::string("checkpoint: ")+message);
}
inline void checkpointEnd(std::istream &in) {
    in>>std::ws;checkpointRequire(in.eof(),"trailing data");
}
inline void checkpointBlob(std::ostream &out,const std::string &s) {
    out<<s.size()<<'\n';out.write(s.data(),s.size());out<<'\n';
    if(!out)throw std::runtime_error("checkpoint section write failed");
}
inline std::string checkpointBlob(std::istream &in) {
    std::uint64_t size=0;checkpointRequire(bool(in>>size),"missing section length");
    checkpointRequire(in.get()=='\n',"invalid section separator");
    const auto here=in.tellg();in.seekg(0,std::ios::end);const auto end=in.tellg();in.seekg(here);
    checkpointRequire(here>=0 && end>=here && size<std::uint64_t(end-here),"truncated section");
    std::string s(std::size_t(size),'\0');in.read(&s[0],size);
    checkpointRequire(bool(in) && in.get()=='\n',"truncated section body");return s;
}
inline std::uint64_t checkpointChecksum(const std::string &s) {
    // Corruption detection, not a cryptographic authenticity certificate.
    std::uint64_t h=14695981039346656037ull;
    for(unsigned char c:s){h^=c;h*=1099511628211ull;}return h;
}
}
#endif
