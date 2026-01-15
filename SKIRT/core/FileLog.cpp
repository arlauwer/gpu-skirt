/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FileLog.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "System.hpp"

////////////////////////////////////////////////////////////////////

FileLog::FileLog() {}

////////////////////////////////////////////////////////////////////

FileLog::FileLog(SimulationItem* parent)
{
    parent->addChild(this);
}

////////////////////////////////////////////////////////////////////

void FileLog::setupSelfBefore()
{
    // Call the setup of the base class first, to ensure the string identifying the process is set.
    Log::setupSelfBefore();

    // Open the log output file
    open();
}

////////////////////////////////////////////////////////////////////

void FileLog::open()
{
    string filepath;
    filepath = find<FilePaths>()->output("log.txt");

    _out = System::ofstream(filepath);
    if (!_out) throw FATALERROR("Could not open the log file " + filepath);
}

////////////////////////////////////////////////////////////////////

namespace
{
    // The strings beginning a message, indexed by level (Info, Warning, Success, Error)
    // NOTE: this depends on the order in the Level enum --> rather dirty
    const char* _messageBegin[] = {"   ", " ! ", " - ", " * "};
}

void FileLog::output(string message, Log::Level level)
{
    if (!_out.is_open() && (level == Level::Warning || level == Level::Error))
    {
        std::unique_lock<std::mutex> lock(_mutex);
        open();
    }

    if (_out.is_open())
    {
        std::unique_lock<std::mutex> lock(_mutex);
        _out << System::timestamp() << _messageBegin[static_cast<size_t>(level)] << message << std::endl;
    }
}

////////////////////////////////////////////////////////////////////
