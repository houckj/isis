#version 1.0
# The initial echo is necessary because the solaris version of sed cannot
# grok input without a trailing newline.
echo `grep "^#define MODULE_[MP]" version.h | sed -e 's/[^0-9]*//' | tr '\012' .` | sed -e 's/.$//'
