#pragma once

#include <exception>
#include <stdexcept>
#include <string>

class UnknownIntegrator : std::runtime_error 
{

public:
  UnknownIntegrator(const std::string& err_msg)
  : std::runtime_error{ "UnknownIntegrator: " + err_msg }
  {}
};

class IntegrationError : std::runtime_error
{
public:
  IntegrationError(const std::string& err_msg)
  : std::runtime_error{ "IntegrationError: " + err_msg }
  {}
};

class EmptyListError : std::runtime_error
{
public:
  EmptyListError(const std::string& err_msg)
  : std::runtime_error{ "EmptyListError: " + err_msg }
  {}
};

class ValueError : std::runtime_error
{
public:
    ValueError(const std::string& err_msg)
    : std::runtime_error{ "ValueError: " + err_msg }
    {}
};
