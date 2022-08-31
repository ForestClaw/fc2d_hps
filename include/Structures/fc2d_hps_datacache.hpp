#ifndef FC2D_HPS_DATACACHE_HPP_
#define FC2D_HPS_DATACACHE_HPP_

#include <mutex>
#include <thread>
#include <map>
#include <string>

#include "Util/GenericSingleton.hpp"

template<class DataType>
class DataCache : public GenericSingleton<DataCache<DataType>> {

public:

    std::map<std::string, DataType> dataMap;

    DataCache() :
        dataMap{}
            {}
    
    DataCache(const std::map<std::string, DataType>& aMap) :
        dataMap(aMap)
            {}

    ~DataCache() {
        dataMap.clear();
    }

};

#endif // FC2D_HPS_DATACACHE_HPP_