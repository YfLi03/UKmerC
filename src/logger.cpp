#include "logger.hpp"
#include <iostream>
#include <numeric>

Logger::Logger() : logstream(new std::ostringstream()) { }


void Logger::flush(std::ostringstream& ss)
{
    flush(ss.str().c_str());
    ss.clear();
    ss.str("");
}


void Logger::flush(char const *label)
{
    std::string slabel(label);
    std::string banner;
    banner.assign(slabel.size(), '=');
    std::cout << "[Logger]" << slabel << "\n    ";
    std::string str = logstream->str();
    std::cout << str << std::endl << std::endl;
    logstream.reset(new std::ostringstream());
}

std::string Logger::readrangestr(size_t pos, size_t count)
{
    std::ostringstream ss;
    ss << "[" << pos << ".." << (pos+count-1) << "] (" << count << " reads)";
    return ss.str();
}
