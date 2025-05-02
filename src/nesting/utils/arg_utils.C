#include "nesting_utils.H"


// format parser for command line options
string_map parse_argv(int argc, char* argv[])
{
  string_map arguments;

  for(int i = 0; i < argc; i++)
  {
    if (argv[i][0] == '-' && argv[i+1][0] != '-')
    {
      arguments[argv[i]] = argv[i+1];
      i++;
    }
  }
  return arguments;
}

// check if the input arguments are valid for the application
int check_argv(const string_map &arguments)
{
  return 0;
}
