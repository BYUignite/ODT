/**
 * @file odtExceptions.h
 * @brief Definitions for ODT-specific error handling
 *
 * Loosely based on Cantera's exceptions and error handling method; see cantera/base/ctexceptions.h
 */

#ifndef ODT_ODTEXCEPTIONS_H
#define ODT_ODTEXCEPTIONS_H

#include "yaml-cpp/yaml.h"
#include "cantera/base/fmt.h"
#include <exception>
#include <string>
#include <sstream>
#include <utility>
#include <cantera/base/ctexceptions.h>

//! Base class for exceptions thrown by ODT classes
class odtError : public std::exception {

public:
    //! Normal Constructor for the odtError base class
    /*!
     * @param procedure Name of the function within which the error was generated. For member functions, this should be written as `ClassName::functionName`. For constructors, this should be `ClassName::ClassName`. Arguments can be specified to disambiguate overloaded functions, such as `ClassName::functionName(int, int)`.
     * @param msg  Descriptive string describing the type of error message. This can be a fmt-style format string (that is, using curly braces to indicate fields), which will be used with additional arguments to generate a formatted error message
     * @param args Arguments which will be used to interpolate the format string
     */
    template <typename... Args>
    odtError(const std::string& procedure, const std::string& msg, const Args&... args) noexcept : procedure_(procedure) {
        if (sizeof...(args) == 0)
            msg_ = msg;
        else {
            msg_ = fmt::format(msg, args...);
        }
    }

    //! Destructor for base class does nothing
    ~odtError() noexcept override = default;

    //! Get a description of the error
    const char* what() const noexcept override;

    //! Method overridden by derived classes to format the error message
    virtual std::string getMessage() const noexcept;

    //! Method overridden by derived classes to indicate their type
    virtual std::string getClass() const noexcept {
        return "odtError";
    }

protected:
    //! Protected default constructor discourages throwing errors containing no information.
    odtError() noexcept = default;

    //! Constructor used by derived classes that override getMessage()
    explicit odtError(const std::string& procedure) noexcept;

    std::string procedure_;                 //! The name of the procedure where the exception occurred
    mutable std::string formattedMessage_;  //!< Formatted message returned by what()

private:
    std::string msg_;                       //!< Message associated with the exception

};

//! An error indicating that the input.yaml file cannot be accessed properly
class inputFileError : public odtError {

public:
    //! Constructor for inputFileError exception
    /*!
     * @param func String name for the function within which the error was
     *             generated.
     */
    explicit inputFileError(const std::string& func, const YAML::BadFile& y) noexcept:
        odtError(func), y_(y) {}

    std::string getMessage() const noexcept override;
    std::string getClass() const noexcept override {
        return "inputFileError";
    }

private:
    YAML::BadFile y_;
};

//! An error indicating that a file cannot be accessed properly
class FileError : public odtError {

public:
    //! Constructor for inputFileError exception
    /*!
     * @param func String name for the function within which the error was
     *             generated.
     */
    explicit FileError(const std::string& func, const std::string& filename) noexcept:
            odtError(func), filename_(filename) {}

    std::string getMessage() const noexcept override;
    std::string getClass() const noexcept override {
        return "FileError";
    }

private:
    std::string filename_;

};

//! An error indicating an invalid input parameter
class inputParamError : public odtError {

public:
    //! Constructor for inputParamError exception
    /*!
     * @@param procedure Name of the function within which the error was generated
     * @param param Name of the invalid parameter
     */
    inputParamError(const std::string& procedure, const std::string& param) noexcept:
        odtError(procedure), param_(param) {}

    std::string getMessage() const noexcept override;
    std::string getClass() const noexcept override {
        return "inputParamError";
    }

private:
    std::string param_;
};

//! An error indicating error thrown by Cantera
class odtCanteraError : public odtError {

public:
    //! Constructor for odtCanteraError exception
    /*!
     * @@param procedure Name of the function within which the error was generated
     * @param param Name of the invalid parameter
     */
    odtCanteraError(const std::string& procedure, const std::string& func, const Cantera::CanteraError &c) noexcept :
        odtError(procedure), func_(func), c_(c) {}

    std::string getMessage() const noexcept override;
    std::string getClass() const noexcept override {
        return "odtCanteraError";
    }

private:
    std::string func_;
    Cantera::CanteraError c_;
};

//! An error indicating error thrown by CVode
class odtCvodeError : public odtError {

public:
    //! Constructor for odtCvodeError exception
    /*!
     * @@param procedure Name of the function within which the error was generated
     * @param param Name of the invalid parameter
     */
    odtCvodeError(std::string& procedure, std::string& param) noexcept :
        odtError(procedure) {}

    std::string getMessage() const noexcept override;
    std::string getClass() const noexcept override {
        return "odtCvodeError";
    }
};

//! Provides a line number
#define XSTR_TRACE_LINE(s) STR_TRACE_LINE(s)

//! Provides a line number
#define STR_TRACE_LINE(s) #s

//! Provides a std::string variable containing the file and line number
#define STR_TRACE (std::string(__FILE__) + ":" + XSTR_TRACE_LINE(__LINE__))

#endif //ODT_ODTEXCEPTIONS_H