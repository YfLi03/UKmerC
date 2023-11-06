#ifndef LOGGER_H_
#define LOGGER_H_

#include <iostream>
#include <memory>
#include <limits>
#include <cstdint>
#include <sstream>

class Logger {
    std::unique_ptr<std::ostringstream> logstream, rootstream;

public:
    Logger();
    void flush(char const *label);
    void flush(std::ostringstream& ss);
    std::ostringstream& operator()() { return *logstream; }
    static std::string readrangestr(size_t pos, size_t count);
};

#endif