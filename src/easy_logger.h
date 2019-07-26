/// @file
/// @brief すぐ使える実用的な簡易ロガー
/// @description LOGI, LOGD, LOGW, LOGE の4つのマクロが定義され、それぞれをログ出力ストリームと見做して LOGI << "hoge"; のように使う。
/// LOGI ... info  レベルのログメッセージ出力（白色）
/// LOGD ... debug レベルのログメッセージ出力（緑色）
/// LOGW ... warn  レベルのログメッセージ出力（黄色）
/// LOGE ... error レベルのログメッセージ出力（赤色）
/// @note
/// (a): USE_UTIL_LOG_EASY_LOGGER_BOOST_CPU_TIMER が定義されている場合、時間計測に boost::timer::cpu_timer を使用する（ boost.timer 関連のライブラリーのリンクが必要 ）
/// (b): USE_UTIL_LOG_EASY_LOGGER_BOOST_CHRONO が定義されている場合、時間計測に boost::chrono::steady_clock を使用する ( boost.chrono 関連のライブラリーのリンクが必要)
/// (c): USE_UTIL_LOG_EASY_LOGGER_STD_CHRONO が定義されている場合、時間計測に std::chrono::steady_clock を使用する（ 外部ライブラリーのリンクは不要だが、処理系によっては分解能が不足する ）
/// (d): (a), (b), (c) の何れも定義されていない場合、時間計測に util::chrono::default_clock を使用する（ 外部ライブラリーは不要、windows処理系でもQueryPerformanceCounterを内部的に使用する ）
/// (f): DISABLE_UTIL_LOG_EASY_LOGGER が定義されている場合、全てのログ出力は事実上無効になる。
/// UTIL_LOG_EASY_LOGGER_GET_{ FUNCTION | FILE | LINE } を事前に定義しておくとユーザー定義の何かそういうものに置き換え可能。（ "" をユーザー定義すれば出力から消す事も可能。 ）
/// UTIL_LOG_EASY_LOGGER_{INFO|DEBUG|WARN|ERROR}_PREFIX を事前に定義しておくと [ info ] 的な部分をユーザー定義に置き換え可能。（ "" をユーザー定義すれば出力から消す事も可能。 ）
//https://qiita.com/usagi/items/d4aec8d3f748f4ba9d6a
//https://github.com/usagi/usagi/blob/master/include/usagi/log/easy_logger.hxx

#pragma once

#ifdef DISABLE_EASY_LOGGER


#define LOGI ::easy_logger::log_null()
#define LOGD ::easy_logger::log_null()
#define LOGW ::easy_logger::log_null()
#define LOGE ::easy_logger::log_null()

#else

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#ifndef EASY_LOGGER_GET_FILE
  #define EASY_LOGGER_GET_FILE __FILENAME__
#endif

#ifndef EASY_LOGGER_GET_LINE
  #define EASY_LOGGER_GET_LINE __LINE__
#endif

#ifndef EASY_LOGGER_GET_FUNCTION
  #if defined( __Clang__ ) || defined( __GNUC__ )
    #define EASY_LOGGER_GET_FUNCTION __FUNCTION__
  #else
    #define EASY_LOGGER_GET_FUNCTION __func__
  #endif
#endif

#ifndef EASY_LOGGER_INFO_PREFIX
  #define EASY_LOGGER_INFO_PREFIX  " [ INFO  ] "
#endif

#ifndef EASY_LOGGER_DEBUG_PREFIX
  #define EASY_LOGGER_DEBUG_PREFIX " [ DEBUG ] "
#endif

#ifndef EASY_LOGGER_WARN_PREFIX
  #define EASY_LOGGER_WARN_PREFIX  " [ WARN ] "
#endif

#ifndef EASY_LOGGER_ERROR_PREFIX
  #define EASY_LOGGER_ERROR_PREFIX " [ ERROR ] "
#endif


#define LOGI ::easy_logger::log_intermediate::make_info \
  ( EASY_LOGGER_GET_FILE \
  , EASY_LOGGER_GET_LINE \
  , EASY_LOGGER_GET_FUNCTION \
  ) << EASY_LOGGER_INFO_PREFIX
#define LOGD ::easy_logger::log_intermediate::make_debug\
  ( EASY_LOGGER_GET_FILE \
  , EASY_LOGGER_GET_LINE \
  , EASY_LOGGER_GET_FUNCTION \
  ) << EASY_LOGGER_DEBUG_PREFIX
#define LOGW ::easy_logger::log_intermediate::make_warn \
  ( EASY_LOGGER_GET_FILE \
  , EASY_LOGGER_GET_LINE \
  , EASY_LOGGER_GET_FUNCTION \
  ) << EASY_LOGGER_WARN_PREFIX
#define LOGE ::easy_logger::log_intermediate::make_error\
  ( EASY_LOGGER_GET_FILE \
  , EASY_LOGGER_GET_LINE \
  , EASY_LOGGER_GET_FUNCTION \
  ) << EASY_LOGGER_ERROR_PREFIX

#endif

#ifndef EASY_LOGGER_OUT
  #define EASY_LOGGER_OUT ::std::cout
#endif

#ifndef EASY_LOGGER_FLUSH
  #define EASY_LOGGER_FLUSH ::std::flush
#endif

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace easy_logger
{
  using namespace std;

  
#ifdef DISABLE_EASY_LOGGER
  
  struct log_null
  {
    template < typename T >
    decltype( auto ) operator<<( const T& ) { return *this; }
  };
  
#else
  
  
  struct log_intermediate
  {
    stringstream buffer;
    
    const char* prefix;
    const char* suffix;
    
    const char* source;
    const size_t line;
    const char* function1;
    
    log_intermediate( log_intermediate&& a )
      : buffer( move( a.buffer ) )
      , source( a.source )
      , line( a.line )
      , function1( a.function1 )
    { }
    
    log_intermediate( const char* source_, const size_t line_, const char* function_ )
      : prefix( p )
      , suffix( s )
      , source( source_ )
      , line( line_ )
      , function1( function_ )
    { }
    
    static log_intermediate make_info ( const char* s, const size_t l, const char* f ) { return log_intermediate( s, l, f ); }
    static log_intermediate make_debug( const char* s, const size_t l, const char* f ) { return log_intermediate( s, l, f ); }
    static log_intermediate make_warn ( const char* s, const size_t l, const char* f ) { return log_intermediate( s, l, f ); }
    static log_intermediate make_error( const char* s, const size_t l, const char* f ) { return log_intermediate( s, l, f ); }
    
    template < typename T >
    stringstream operator<<( const T& in )
    {
      return buffer << in;
    }
    
    ~log_intermediate() noexcept
    {
      //stringstream s;
      //s << std::left << std::setw(30) << function << " " << line
      //  << prefix
      //  << buffer.str()
      //  << " "
      //  << suffix
      //  << flush;
      //s << source << '\t' << line << '\t' << function
      //  << prefix
      //  << buffer.str()
      //  << "\t"
      //  << suffix
      //  << endl
      //  ;
      //EASY_LOGGER_OUT << s.str();
      std::cout << std::left << std::setw(30) << function1 << " " << line << prefix << buffer.str() << " " << suffix << endl;
      //EASY_LOGGER_OUT << s.str() << EASY_LOGGER_FLUSH;
    }
  };

#endif

}
