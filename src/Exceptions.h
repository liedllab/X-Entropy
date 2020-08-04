#pragma once

#include <exception>
#include <string>

class UnknownIntegrator : std::exception 
{
private:
  std::string m_error_message;

public:
  UnknownIntegrator(const std::string& err_msg)
  : m_error_message{ "Error: Unknown Integrator: " + err_msg }
  {

  }

  virtual const char *what() const throw() override
  {
    return m_error_message.c_str();
  }
};

class IntegrationError : std::exception
{
private:
  std::string m_error_message;
public:
  IntegrationError(const std::string& err_msg)
  : m_error_message{ "Integration Error: " + err_msg }
  {

  }

  virtual const char *what() const throw() override
  {
    return m_error_message.c_str();
  }
};