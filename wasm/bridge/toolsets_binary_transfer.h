#pragma once
#include <nlohmann/json.hpp>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

// ZST1: little-endian tagged tree. Numeric leaves are contiguous int32/float64
// blocks; the public object structure is independent of JSON text formatting.
namespace zspaceTransfer {
using Json=nlohmann::json;
constexpr uint32_t magic=0x3154535a;
constexpr size_t limit=256u*1024u*1024u;
struct Reader {
    const uint8_t* data;size_t size,offset=0;
    template<class T>T scalar(){if(sizeof(T)>size-offset)throw std::invalid_argument("Truncated toolsets buffer.");T v;std::memcpy(&v,data+offset,sizeof(T));offset+=sizeof(T);return v;}
    std::string string(){auto n=scalar<uint32_t>();if(n>size-offset)throw std::invalid_argument("Truncated toolsets string.");std::string s(reinterpret_cast<const char*>(data+offset),n);offset+=n;return s;}
    Json value(unsigned depth=0){
        if(depth>64)throw std::invalid_argument("Toolsets nesting limit exceeded.");
        auto tag=scalar<uint8_t>();
        switch(tag){
        case 0:return nullptr;case 1:return false;case 2:return true;
        case 3:{double n=scalar<double>();if(!std::isfinite(n)||std::abs(n)>9007199254740991.0||std::floor(n)!=n)throw std::invalid_argument("Invalid toolsets integer.");return int64_t(n);}
        case 4:return scalar<double>();case 5:return string();
        case 6:{auto n=scalar<uint32_t>();if(n>(size-offset)/5)throw std::invalid_argument("Invalid toolsets object size.");Json out=Json::object();for(uint32_t i=0;i<n;i++){auto key=string();if(out.contains(key))throw std::invalid_argument("Duplicate toolsets key.");out[key]=value(depth+1);}return out;}
        case 7:case 8:case 9:{auto n=scalar<uint32_t>();size_t width=tag==8?4:tag==9?8:1;if(n>(size-offset)/width)throw std::invalid_argument("Invalid toolsets array size.");Json out=Json::array();auto& a=out.get_ref<Json::array_t&>();a.reserve(n);for(uint32_t i=0;i<n;i++)a.push_back(tag==8?Json(scalar<int32_t>()):tag==9?Json(scalar<double>()):value(depth+1));return out;}
        default:throw std::invalid_argument("Unknown toolsets buffer tag.");
        }
    }
};
inline Json decode(const void* data,size_t size){
    if(!data||size<5||size>limit)throw std::invalid_argument("Invalid toolsets buffer size.");
    Reader r{static_cast<const uint8_t*>(data),size};if(r.scalar<uint32_t>()!=magic)throw std::invalid_argument("Unsupported toolsets transfer version.");
    auto out=r.value();if(r.offset!=size)throw std::invalid_argument("Trailing toolsets buffer bytes.");return out;
}
struct Writer {
    std::vector<uint8_t> bytes;
    template<class T>void scalar(T v){if(bytes.size()>limit-sizeof(T))throw std::invalid_argument("Toolsets result exceeds buffer limit.");auto p=reinterpret_cast<const uint8_t*>(&v);bytes.insert(bytes.end(),p,p+sizeof(T));}
    void string(const std::string& s){if(s.size()>limit-4||bytes.size()>limit-4-s.size())throw std::invalid_argument("Toolsets string exceeds limit.");scalar<uint32_t>(uint32_t(s.size()));bytes.insert(bytes.end(),s.begin(),s.end());}
    void value(const Json& v,unsigned depth=0){
        if(depth>64)throw std::invalid_argument("Toolsets nesting limit exceeded.");
        if(v.is_null()){scalar<uint8_t>(0);return;}
        if(v.is_boolean()){scalar<uint8_t>(v.get<bool>()?2:1);return;}
        if(v.is_number()){scalar<uint8_t>(v.is_number_integer()?3:4);scalar<double>(v.get<double>());return;}
        if(v.is_string()){scalar<uint8_t>(5);string(v.get_ref<const std::string&>());return;}
        if(v.is_object()){scalar<uint8_t>(6);scalar<uint32_t>(uint32_t(v.size()));for(auto it=v.begin();it!=v.end();++it){string(it.key());value(it.value(),depth+1);}return;}
        if(v.is_array()){
            bool numeric=true,integer=true;for(const auto& x:v){if(!x.is_number()){numeric=false;break;}auto n=x.get<double>();if(!x.is_number_integer()||n<INT32_MIN||n>INT32_MAX)integer=false;}
            auto tag=numeric?(integer?8:9):7;scalar<uint8_t>(uint8_t(tag));scalar<uint32_t>(uint32_t(v.size()));
            for(const auto& x:v){if(tag==8)scalar<int32_t>(x.get<int32_t>());else if(tag==9)scalar<double>(x.get<double>());else value(x,depth+1);}return;
        }
        throw std::invalid_argument("Unsupported toolsets value.");
    }
};
inline std::vector<uint8_t> encode(const Json& value){Writer w;w.scalar<uint32_t>(magic);w.value(value);return std::move(w.bytes);}
}
