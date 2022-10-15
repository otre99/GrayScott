#ifndef SIMPLE_CMD_PARSER_H
#define SIMPLE_CMD_PARSER_H

#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

class SimpleCmdParser {
 public:
  SimpleCmdParser();
  bool Parse(int argc, char *argv[], const std::string &program_info = {});
  void SetDefaultValue(const std::string &key, const std::string &value,
                       const std::string &info = {});
  bool Has(const std::string &key);
  void PrintHelp();

  template <typename dtype>
  dtype Get(const std::string &key);

 private:
  std::map<std::string, std::string> cmds_;
  std::map<std::string, std::string> help_info_;
  std::string program_info_;
};

SimpleCmdParser::SimpleCmdParser() { cmds_.clear(); }

bool SimpleCmdParser::Parse(int argc, char *argv[], const std::string &info) {
  program_info_ = info;
  for (int i = 1; i < argc; i += 2) {
    try {
      const auto key = std::string(argv[i]);
      if (!this->Has(key)) {
        std::cerr << "Unknow parameter '" << key << "'" << std::endl;
        return false;
      }
      cmds_[key] = std::string(argv[i + 1]);
    } catch (...) {
      return false;
    }
  }
  return true;
}

void SimpleCmdParser::SetDefaultValue(const std::string &key,
                                      const std::string &value,
                                      const std::string &info) {
  cmds_[key] = value;
  help_info_[key] = info;
}

void SimpleCmdParser::PrintHelp() {
  std::cout << "Info: " << program_info_ << std::endl;
  std::cout << "Arguments:" << std::endl;
  for (const auto &o : cmds_) {
    std::cout << std::setw(12) << o.first << "    Help: " << help_info_[o.first]
              << ". Default ( " << o.second << " )" << std::endl;
  }
}

bool SimpleCmdParser::Has(const std::string &key) {
  return cmds_.find(key) != cmds_.end();
}

template <typename dtype>
dtype SimpleCmdParser::Get(const std::string &key) {
  dtype r = {};
  if (!Has(key)) {
    std::cerr << "Error. Not found any the argument " << key << std::endl;
    return r;
  }
  std::istringstream in(cmds_[key]);
  in >> r;
  return r;
}

#endif  // SIMPLE_CMD_PARSER_H
