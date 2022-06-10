/**
 * @file odtExceptions.cc
 * @brief Definitions for ODT-specific error handling
 *
 * Loosely based on Cantera's exceptions and error handling method; see cantera/base/ctexceptions.cpp
 */

#include "odtExceptions.h"

// *** Exceptions ***

static const char* line = "-----------------------------------------------------------------------\n";

odtError::odtError(const std::string& procedure) noexcept : procedure_(procedure)
{
}

const char* odtError::what() const noexcept
{
    try {
        formattedMessage_ = "\n";
        formattedMessage_ += line;
        formattedMessage_ += getClass();
        if (procedure_.size()) {
            formattedMessage_ += " thrown by " + procedure_;
        }
        formattedMessage_ += "\n" + getMessage();
        if (formattedMessage_.compare(formattedMessage_.size()-1, 1, "\n")) {
            formattedMessage_.append("\n");
        }
        formattedMessage_ += line;
    } catch (...) {
        // Something went terribly wrong and we couldn't even format the message.
    }
    return formattedMessage_.c_str();
}

std::string odtError::getMessage() const noexcept
{
    return msg_;
}

std::string inputFileError::getMessage() const noexcept
{
    return fmt::format("Unable to access YAML input file.");
}

std::string FileError::getMessage() const noexcept
{
    return fmt::format("File access error: {}", filename_);
}

std::string inputParamError::getMessage() const noexcept
{
    return fmt::format("Invalid input parameter: {}", param_);
}

std::string odtCanteraError::getMessage() const noexcept
{
    return fmt::format("CanteraError thrown from {} {}", func_, c_.what());
}

std::string odtCvodeError::getMessage() const noexcept
{
    return fmt::format("Error thrown by CVode");
}