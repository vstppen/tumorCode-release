#include <cstdlib>
#include <string>

inline bool clear_directory(const std::string& path) {
    std::string command = "rm -rf " + path + "/*";  // Linux/Mac
    return system(command.c_str()) == 0;
}