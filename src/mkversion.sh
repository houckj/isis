grep "^#define MODULE_[MP]" version.h | sed -e 's/[^0-9]*//' | tr '\012' . | sed -e 's/.$//'
