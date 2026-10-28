// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/setup_registry.hpp
//
// Last modified : 10/2026
//
// Copyright(C) by :
// - Geoffroy Lesur <geoffroy.lesur@univ-grenoble-alpes.fr> (IPAG/UGA/CNRS - 2021 - 2026)
// and other code contributors
//
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef SETUP_REGISTRY_HPP_
#define SETUP_REGISTRY_HPP_

#include <functional>
#include <memory>
#include <utility>  // for std::move
#include "setup.hpp"

class SetupRegistry {
 public:
    using SetupFactory = std::function<std::unique_ptr<SetupBase>(
        Input&, Grid&, DataBlock&, Output&)>;

    // Register the user's setup factory (only one allowed per binary)
    static void Register(SetupFactory factory) {
        auto& instance = GetInstance();
        if (instance.factory) {
            IDEFIX_ERROR("Multiple Setup registrations detected. "
                         "Only one REGISTER_SETUP() call is allowed per binary.");
        }
        instance.factory = std::move(factory);
    }

    // Create the setup: uses registered factory or falls back to DefaultSetup
    static std::unique_ptr<SetupBase> Create(
        Input& input, Grid& grid, DataBlock& data, Output& output) {
        auto& instance = GetInstance();
        if (instance.factory) {
            return instance.factory(input, grid, data, output);
        }
        // Fallback to default (backward compatible with old Setup class)
        return std::make_unique<DefaultSetup>(input, grid, data, output);
    }

 private:
    SetupFactory factory;

    static SetupRegistry& GetInstance() {
        static SetupRegistry instance;  // Meyer's singleton
        return instance;
    }

    SetupRegistry() = default;
    ~SetupRegistry() = default;
    SetupRegistry(const SetupRegistry&) = delete;
    SetupRegistry& operator=(const SetupRegistry&) = delete;
};

// Macro for one-line registration in user's setup.cpp
#define REGISTER_SETUP(SetupClass) \
    namespace { \
        struct SetupRegistrar { \
            SetupRegistrar() { \
                SetupRegistry::Register( \
                    [](Input& input, Grid& grid, DataBlock& data, Output& output) { \
                        return std::make_unique<SetupClass>( \
                            input, grid, data, output); \
                    }); \
            } \
        } setup_registrar; \
    }  // namespace

#endif  // SETUP_REGISTRY_HPP_
