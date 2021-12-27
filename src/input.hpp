// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef INPUT_HPP_
#define INPUT_HPP_

#include <string>
#include <map>
#include <vector>

#include "idefix.hpp"

using IdefixParamContainer = std::vector<std::string>;
using IdefixBlockContainer = std::map<std::string,IdefixParamContainer>;
using IdefixInputContainer = std::map<std::string,IdefixBlockContainer>;

class Input {
 public:
  // Constructor from a file
  Input (int, char ** );
  void PrintParameters();

  // Accessor to input parameters
  // the parameters are always: BlockName, EntryName, ParameterNumber (starting from 0)
  std::string GetString(std::string, std::string, int); ///< Read a string from the input file
  real GetReal(std::string, std::string, int);          ///< Read a real number from the input file
  int GetInt(std::string, std::string, int);            ///< Read an integer from the input file
  int CheckEntry(std::string, std::string);             ///< Check that a block+entry is present
                                                        ///< in the input file
  bool CheckBlock(std::string);                         ///< check that whether a block is defined
                                                        ///< in the input file
  bool CheckForAbort();                                 // have we been asked for an abort?
  void CheckForStopFile();                              // have we been asked for an abort from
                                                        // a stop file?

  Input();
  void PrintLogo();

  bool restartRequested{false};       //< Should we restart?
  int  restartFileNumber;             //< if yes, from which file?

  static bool abortRequested;         //< Did we receive an abort signal (USR2) from the system?

  bool tuningRequested{false};        //< whether the user has asked for loop-tuning

  int maxCycles{-1};                   //< whether we should perform a maximum number of cycles

  bool forceNoWrite{false};           //< explicitely disable all writes to disk

 private:
  std::string inputFileName;
  IdefixInputContainer  inputParameters;
  void ParseCommandLine(int , char **argv);
  static void signalHandler(int);
  std::vector<std::string> getDirectoryFiles();
  std::string getFileExtension(const std::string file_name);
};

#endif // INPUT_HPP_
