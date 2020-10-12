#pragma once

#include <exception>
#include <string>

class UnknownIntegrator : std::runtime_error 
{

public:
  UnknownIntegrator(const std::string& err_msg)
  : std::runtime_error{ "Error: Unknown Integrator: " + err_msg }
  {}
};

class IntegrationError : std::runtime_error
{
public:
  IntegrationError(const std::string& err_msg)
  : std::runtime_error{ "Integration Error: " + err_msg }
  {}
};

class EmptyListError : std::runtime_error
{
public:
  EmptyListError(const std::string& err_msg)
  : std::runtime_error{ "Empty List Error: " + err_msg }
  {}
};
