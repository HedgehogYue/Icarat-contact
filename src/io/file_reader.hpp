///  @file  file_reader.hpp
///  @author	Daiki Watanabe
///  @date		January 22, 2021.
/// Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
/// This software is released under the MIT License.
/// http:opensource.org/licenses/mit-license.php

#pragma once
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace icarat
{
//  read text and convert to vector of strings
template <class T> std::vector<std::string> read_file(T filePass)
{
  std::string line;
  std::vector<std::string> lines;

  std::ifstream input_file(filePass);
  if(!input_file.is_open())
  {
    std::cerr << "Could not open the file - '" << filePass << "'" << std::endl;
    exit(1);
  }

  while(getline(input_file, line))
    lines.push_back(line);

  input_file.close();

  return lines;
}
} // namespace icarat