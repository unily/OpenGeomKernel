#ifndef EXTENSION_EXPORT_H
#define EXTENSION_EXPORT_H

#ifdef _WIN32
#ifdef  Extension_EXPORTS
#define EXTENSION_EXPORT __declspec(dllexport)
#else
#define EXTENSION_EXPORT __declspec(dllimport)
#endif
#else
#define EXTENSION_EXPORT __attribute__((visibility("default")))
#endif

#endif

