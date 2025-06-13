#pragma once

using namespace AMDiS;
using namespace Dune::Functions::BasisFactory;
using namespace Dune::Indices;

template <typename Writer>
void writePVDFile (std::vector<double> const& timesteps,
                   std::string const& name,
                   std::string const& path,
                   Writer& vtkWriter,
                   int gridDim)
{
    /* remember current time step */
    unsigned int count = timesteps.size();

    std::ofstream pvdFile;
    pvdFile.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                       std::ios_base::eofbit);
    std::string pvdname = path + "/" + name + ".pvd";
    pvdFile.open(pvdname.c_str());
    pvdFile << "<?xml version=\"1.0\"?> \n"
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\""
            << Dune::VTK::getEndiannessString() << "\"> \n"
            << "<Collection> \n";
    double timestepSize = 10;
    Parameters::get("adapt->timestep",timestepSize);

    // filename
    std::string fullname;
    for (unsigned int i=0; i<count; i++)
    {
        std::stringstream n;
        n.fill('0');
        n << name << "-" << std::setw(5) << i;

        fullname = n.str();

        pvdFile << "<DataSet timestep=\"" << timesteps[i]
                << "\" group=\"\" part=\"0\" name=\"\" file=\"_piecefiles/"
                << fullname << (gridDim == 1 ? ".vtp" : ".vtu") << "\"/> \n";
    }
    pvdFile << "</Collection> \n"
            << "</VTKFile> \n" << std::flush;
    pvdFile.close();

    //write the file of the actual time step
    vtkWriter.write(path + "/_piecefiles/" + fullname);
}


template <class Writer, class T>
void addWriterData(Writer& writer, T const& arg, std::string name, int dim) {
    if (dim == 1) {
        writer.addVertexData(arg, Dune::VTK::FieldInfo{name, Dune::VTK::FieldInfo::Type::scalar, 1});
    }
    if (dim == 2) {
        writer.addVertexData(arg, Dune::VTK::FieldInfo{name, Dune::VTK::FieldInfo::Type::vector, 2});
    }
}
template <class Writer, class T, class... Types>
void addWriterData(Writer& writer, T const& firstArg, std::string firstName, int firstDim, Types... args) {
    if (firstDim == 1) {
        writer.addVertexData(firstArg, Dune::VTK::FieldInfo{firstName, Dune::VTK::FieldInfo::Type::scalar, 1});
    }
    if (firstDim == 2) {
        writer.addVertexData(firstArg, Dune::VTK::FieldInfo{firstName, Dune::VTK::FieldInfo::Type::vector, 2});
    }
    addWriterData(writer, args...);
}

template <class P, class BW, class KV, class DV,class ExtensionManager>
void writeFiles(P& nschProb,
                BW& bulkWriter,
                KV const& kappaVecOld,
                DV const& y0,
                std::vector<double>& timesteps,
                AdaptInfo& adaptInfo,
                std::string const& path,
                ExtensionManager& extensionManager = {},
                int initial = 0) {
    Dune::Timer t;

    auto const& surfacePreBasis = nschProb.globalBasis()->preBasis().subPreBasis(_phiS);
    auto const& surfaceGrid = surfacePreBasis.surfaceGrid();

    Dune::VTKWriter surfWriter(surfaceGrid->leafGridView(), Dune::VTK::conforming);
    auto mu = nschProb.solution(_mu);
    auto phi = nschProb.solution(_phi);
    auto p = nschProb.solution(_p);
    auto v = nschProb.solution(_v);

    if (!initial) {
        auto axi = Parameters::get<int>("axisymmetric").value_or(0);
        auto Ks = Parameters::get<double>("areaShear").value_or(0);

        // update distances, lambda2 and normal for writeFiles
        updateDistances(nschProb);
        if (axi && Ks) updateLambda2(nschProb, y0);
        nschProb.updateNormal();
    }

    // store the desired solutions in the writer
    if (initial) {
        if (adaptInfo.time() > 0) {
            nschProb.globalBasis()->update(nschProb.gridView());
            nschProb.updateSurfaceBasis();
        }
        addWriterData(bulkWriter,
                      phi, "phi", 1,
                      p, "p", 1,
                      v, "v", 2,
                      mu, "mu", 1,
                      valueOf(*nschProb.distancesDOF()), "lambda1Old", 1,
                      valueOf(*nschProb.lambda2()), "lambda2Old", 1,
                      valueOf(*nschProb.deltaP()), "deltaP", 1,
                      valueOf(*nschProb.deltaPhi()), "deltaPhi", 1,
                      nschProb.nuPhase(), "nuPhase", 1,
                      valueOf(*nschProb.normalDOF()), "normal", 2);
    }

    // store surface normals and surface curvature into surface writer
    auto kOld = surfaceGridFct2D(kappaVecOld, surfacePreBasis);
    DOFVector normalDOFS(nschProb.gridView(), power<2>(Dune::Surfacebasis::surfaceLagrange<1>(nschProb.partitionsRef()),
                                                       flatInterleaved()));
    valueOf(normalDOFS) << valueOf(*nschProb.normalDOF());
    auto n = surfaceGridFct2D(valueOf(normalDOFS), surfacePreBasis);

    DOFVector muS(nschProb.gridView(), Dune::Surfacebasis::surfaceLagrange<1>(nschProb.partitionsRef()));
    surfaceValueOf(muS,nschProb.partitionsRef(),0) << nschProb.solution(_mu);
    auto mus = surfaceGridFct(valueOf(muS), surfacePreBasis);

    addWriterData(surfWriter,
                  kOld,"KappaTimesNOld",2,
                  n,"normal",2);

    auto const &spb = nschProb.globalBasis()->preBasis().subPreBasis(_phiS);

    auto coordsSurface = surfaceGridFct2D(nschProb.solution(_xS), spb,_xS);
    auto laplaceSurface = surfaceGridFct2D(nschProb.solution(_kappaVec), spb,_kappaVec);
    auto lambda1 = surfaceGridFct(nschProb.solution(_lambda1), spb,_lambda1);
    auto lambda2 = surfaceGridFct(nschProb.solution(_lambda2), spb,_lambda2);
    auto force = surfaceGridFct2D(nschProb.solution(_fShell), spb,_fShell);
    auto vS = surfaceGridFct2D(nschProb.solution(_vS), spb,_vS);
    auto phiS = surfaceGridFct(nschProb.solution(_phiS), spb,_phiS);

    // add data to the writer (source functions must still exist when write() function is called later!)
    addWriterData(surfWriter, coordsSurface,"xS",2,
                              laplaceSurface,"KappaTimesN",2,
                              lambda1,"lambda1",1,
                              lambda2,"lambda2",1,
                              force,"fShell",2,
                              vS,"vS",2,
                              phiS,"phiS",1,
                              mus,"muS",1);

    using SurfaceGridFct = decltype(phiS);
    using SurfaceGridFct2D = decltype(vS);
    // write files for extensions
    extensionManager.writeFiles(); //create the output data GridFunctions
    auto outputSurfaceFields = extensionManager.getOutputSurfaceFields();
    // iterate over all outputFields
    for (int i = 0; i < outputSurfaceFields.size(); ++i) {
        // iterate over a specific extension
        for (const auto& [name, value] : outputSurfaceFields[i]) {
            value->addToSurfaceWriter(surfWriter, name);
        }
    }

    if (initial) {
        auto outputBulkFields = extensionManager.getOutputBulkFields();
        // iterate over all outputFields
        for (int i = 0; i < outputBulkFields.size(); ++i) {
            // iterate over a specific extension
            for (const auto &[name, value]: outputBulkFields[i]) {
                value->addToWriter(bulkWriter, name);
            }
        }
    }

    // write the files of the 0-th step of the simulation to .pvd files in th output folder
    if (initial && adaptInfo.time() == 0 || !initial) {
        // write (full) solution to (.pvd) file
        timesteps.push_back(adaptInfo.time());
        writePVDFile<decltype(bulkWriter)>(timesteps, "bulk",  path, bulkWriter, 2);
        writePVDFile<decltype(surfWriter)>(timesteps, "surface",  path, surfWriter, 1);
    } else {
        //write the file of the actual time step
        bulkWriter.write( path + "/bulkAfterLastRemeshing");
        surfWriter.write( path + "/surfaceAfterLastRemeshing");
    }

    AMDiS::info(2, "write File step needed {} seconds", t.elapsed());
}
