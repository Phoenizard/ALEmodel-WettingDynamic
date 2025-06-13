//
// Created by mmokbel on 15.05.25.
//
#ifndef ALEMODELAMDIS2_EXTENSIONINTERFACE_HPP
#define ALEMODELAMDIS2_EXTENSIONINTERFACE_HPP

using Writer = typename Dune::VTKWriter<NavierStokesCahnHilliard<CurvedGrid,PBF>::GridView>;
static const int dim = GridView::dimension - 1; // surface grid dimension
using WriterSurface = typename Dune::VTKWriter<Dune::FoamGrid<dim,GridView::dimension>::LeafGridView>;

/// an interface for output fields for the file writer to work properly
class AbstractOutputField {
public:
    virtual ~AbstractOutputField() = default;
    virtual void addToWriter(Writer& writer, const std::string& name) const = 0;
    virtual void addToSurfaceWriter(WriterSurface& writer, const std::string& name) const = 0;
};

/// implementation of AbstractOutputField for bulk output
template <class DOFVectorType, int compDim>
class OutputField : public AbstractOutputField {
public:
    explicit OutputField(std::shared_ptr<DOFVectorType> vec) : vec_(std::move(vec)) {}

    void addToWriter(Writer& writer, const std::string& name) const override {
        addWriterData(writer, valueOf(*vec_), name, compDim);
    }

    void addToSurfaceWriter(WriterSurface& writer, const std::string& name) const override {}

private:
    std::shared_ptr<DOFVectorType> vec_;
};

/// implementation of AbstractOutputField for surface output
template <class DOFVectorType, int compDim>
class OutputSurfaceField : public AbstractOutputField {
public:
    explicit OutputSurfaceField(std::shared_ptr<DOFVectorType> vec) : vec_(std::move(vec)) {}

    void addToSurfaceWriter(WriterSurface& writer, const std::string& name) const override {
        addWriterData(writer, *vec_, name, compDim);
    }

    void addToWriter(Writer& writer, const std::string& name) const override {}

private:
    std::shared_ptr<DOFVectorType> vec_;
};

/// the interface to implement extensions of the membrane wetting model
class ExtensionInterface {
public:
    /// method for initial creation of data structures needed throughout the simulation
    virtual void initData(AdaptInfo &adaptInfo) = 0;

    /// method that is called before the system is assembled and solved,
    /// fill with updates of data structures, variables etc.
    virtual void initTimestep(AdaptInfo &adaptInfo) = 0;

    /// method that is called after the system has been solved,
    /// fill with extra output calculations, updates, etc.
    virtual void closeTimestep(AdaptInfo &adaptInfo) = 0;

    /// method that is called once at the beginning of the simulation to add
    /// your own operators to the NSCHProb
    virtual void fillOperators(AdaptInfo &adaptInfo) = 0;

    /// method that is called once at the beginning of the simulation to add
    /// your own boundary conditions to the NSCHProb
    virtual void fillBoundaryConditions(AdaptInfo &adaptInfo) = 0;

    /** method that is called during writeFiles() in each time step
        * you need to fill outputSurfaceFields_ and/or outputBulkFields_ with GridFunctions
        * that you would like to be included in the output files surface.pvd and bulk.pvd
        * these will then be automatically collected and written
        * for bulk.pvd, therefore, create instances of types OutputSurfaceField or OutputField
        * by putting your gridFunction into a std::shared_ptr<OutputField> or
        * std::shared_ptr<OutputSurfaceField>
    **/
    virtual void writeFiles() = 0;

    virtual ~ExtensionInterface() = default;

    auto getOutputSurfaceFields() {
        return outputSurfaceFields_;
    }

    auto getOutputBulkFields() {
        return outputBulkFields_;
    }

protected:
    std::unordered_map<std::string, std::shared_ptr<AbstractOutputField>> outputSurfaceFields_;
    std::unordered_map<std::string, std::shared_ptr<AbstractOutputField>> outputBulkFields_;
};

/// registry for automatic inclusion of all existing extensions
template <typename NSCHProb>
class ExtensionRegistry {
public:
    using Factory = std::function<std::shared_ptr<ExtensionInterface>(NSCHProb&)>;

    static ExtensionRegistry& instance() {
        static ExtensionRegistry registry;
        return registry;
    }

    void registerFactory(std::string const& name, Factory factory) {
        factories_[name] = std::move(factory);
    }

    std::map<std::string, Factory> const& getFactories() const {
        return factories_;
    }

private:
    std::map<std::string, Factory> factories_;
};

/// helper for automatic registry
template <typename NSCHProb, typename ExtensionType>
struct ExtensionAutoRegister {
    explicit ExtensionAutoRegister(std::string const& name) {
        ExtensionRegistry<NSCHProb>::instance().registerFactory(name, [](NSCHProb& prob) {
            return std::make_shared<ExtensionType>(prob);
        });
    }
};

/// macro for registration
#define REGISTER_EXTENSION(NSCHProbType, Name, Type) \
    static ExtensionAutoRegister<NSCHProbType, Type<NSCHProbType>> reg_##Type(Name)



#endif //ALEMODELAMDIS2_EXTENSIONINTERFACE_HPP