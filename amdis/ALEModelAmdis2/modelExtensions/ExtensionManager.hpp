//
// Created by mmokbel on 15.05.25.
//

#ifndef ALEMODELAMDIS2_EXTENSIONMANAGER_HPP
#define ALEMODELAMDIS2_EXTENSIONMANAGER_HPP


template <class NSCHProb>
class ExtensionManager {
    using ExtensionPtr = std::shared_ptr<ExtensionInterface>;
    std::vector<ExtensionPtr> extensions_;
    std::unordered_map<std::string, ExtensionPtr> nameMap;

public:
    explicit ExtensionManager(NSCHProb& prob) {
        auto const &factories = ExtensionRegistry<NSCHProb>::instance().getFactories();
        for (auto const &[name, factory]: factories) {
            if (Parameters::get<int>("extensions->" + name).value_or(0)) {
                addExtension(name,factory(prob));
            }
        }
    }

    void addExtension(const std::string& name, ExtensionPtr ext) {
        nameMap[name] = ext;
        extensions_.push_back(ext);
    }

    ExtensionPtr operator[](size_t i) {
        return extensions_.at(i); // with bounds-check
    }

    ExtensionPtr operator[](const std::string& name) {
        auto it = nameMap.find(name);
        if (it == nameMap.end())
            throw std::out_of_range("Extension not found: " + name);
        return it->second;
    }

    size_t size() const {
        return extensions_.size();
    }

    bool contains(const std::string& name) {
        auto it = nameMap.find(name);
        if (it == nameMap.end())
            return false;
        return true;
    }

    void initData(AdaptInfo& adaptInfo) {
        for (auto& extension : extensions_)
            extension->initData(adaptInfo);
    }

    void initTimestep(AdaptInfo& adaptInfo) {
        for (auto& extension : extensions_)
            extension->initTimestep(adaptInfo);
    }

    void closeTimestep(AdaptInfo& adaptInfo) {
        for (auto& extension : extensions_)
            extension->closeTimestep(adaptInfo);
    }

    void fillOperators(AdaptInfo& adaptInfo) {
        for (auto& extension : extensions_)
            extension->fillOperators(adaptInfo);
    }

    void fillBoundaryConditions(AdaptInfo& adaptInfo) {
        for (auto& extension : extensions_)
            extension->fillBoundaryConditions(adaptInfo);
    }

    void writeFiles() {
        for (auto& extension : extensions_)
            extension->writeFiles();
    }

    auto getOutputSurfaceFields() {
        std::vector<std::unordered_map<std::string, std::shared_ptr<AbstractOutputField>>> outputFieldsList;
        for (auto& extension : extensions_)
            outputFieldsList.push_back(extension->getOutputSurfaceFields());

        return outputFieldsList;
    }

    auto getOutputBulkFields() {
        std::vector<std::unordered_map<std::string, std::shared_ptr<AbstractOutputField>>> outputFieldsList;
        for (auto& extension : extensions_)
            outputFieldsList.push_back(extension->getOutputBulkFields());

        return outputFieldsList;
    }
};

#endif //ALEMODELAMDIS2_EXTENSIONMANAGER_HPP
