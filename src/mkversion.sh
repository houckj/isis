grep "^#define MODULE_[MP]" version.h | sed -e 's/[^0-9]*//' | tr '\n' . | sed -e 's/.$//'
