#ifndef CRAMS_UTILS_LOGGING_H_
#define CRAMS_UTILS_LOGGING_H_

#include <plog/Appenders/ConsoleAppender.h>
#include <plog/Appenders/RollingFileAppender.h>
#include <plog/Formatters/CsvFormatter.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Init.h>
#include <plog/Log.h>

#include "crams/utils/git_revision.h"

inline void log_startup_information(bool quiet = false) {
  if (plog::get() != nullptr) return;  // already initialized

  if (quiet) {
    // Severity none: all LOG* macros are no-ops; no file is opened.
    static plog::ConsoleAppender<plog::TxtFormatter> nullAppender;
    plog::init(plog::none, &nullAppender);
    return;
  }

  // Unlimited rolling: no rotation for scientific runs where full history matters.
  // Log is written to output/cramslog.csv relative to the working directory.
  static plog::RollingFileAppender<plog::CsvFormatter> fileAppender("output/cramslog.csv");
  static plog::ConsoleAppender<plog::TxtFormatter> consoleAppender;

#ifdef DEBUG
  plog::init(plog::debug, &fileAppender).addAppender(&consoleAppender);
#else
  plog::init(plog::info, &fileAppender).addAppender(&consoleAppender);
#endif

  LOGI << "Welcome to CRAMS version " << CRAMS::get_version();
  LOGI << "built on " << __DATE__ << " at " << __TIME__;
  LOGI << "git SHA1: " << CRAMS::git_sha1();
  if (CRAMS::git_has_local_changes())
    LOGW << "working tree has uncommitted changes — results may not be reproducible";
  else
    LOGI << "working tree is clean";
}

#endif  // CRAMS_UTILS_LOGGING_H_
