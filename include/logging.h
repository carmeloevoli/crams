#ifndef INCLUDE_LOGGING_H_
#define INCLUDE_LOGGING_H_

#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Appenders/ConsoleAppender.h>
#include <plog/Appenders/RollingFileAppender.h>
#include <plog/Formatters/CsvFormatter.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Init.h>
#include <plog/Log.h>

void log_startup_information() {
  static plog::RollingFileAppender<plog::CsvFormatter> fileAppender("output/cramslog.csv", 8000, 3);
  static plog::ConsoleAppender<plog::TxtFormatter> consoleAppender;
#ifdef DEBUG
  plog::init(plog::debug, &fileAppender).addAppender(&consoleAppender);
#else
  plog::init(plog::info, &fileAppender).addAppender(&consoleAppender);
#endif
  //   static plog::ConsoleAppender<plog::TxtFormatter> consoleAppender;
  // #ifdef DEBUG
  //   plog::init(plog::debug, &consoleAppender);
  // #else
  //   plog::init(plog::info, &consoleAppender);
  // #endif
  LOGI << "Welcome to CRAMS version " << get_version();
  LOGI << "was built on " << __DATE__ << " at " << __TIME__;
  LOGI << "git version is " << git_sha1();
  LOGW << "has local changes " << git_has_local_changes();
}

#endif  // INCLUDE_LOGGING_H_