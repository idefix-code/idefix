// ***********************************************************************************
// Idefix MHD astrophysical code
//
// Source file src/setup_registry.hpp
//
// Last modified : 11/2026
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
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "setup.hpp"

// Registry of user-defined setups (OOP style, inheriting from SetupBase).
//
// Several setups can be registered in a single binary, each one under a
// distinct name. The setup actually used is selected at run time with the
// [Setup] name entry of the input file:
//
//  - name = <registered name> : instantiate that setup
//  - name = legacy            : use the legacy Setup class (weak symbols)
//  - no name                  : use the legacy Setup class (with a warning
//                               when at least one setup is registered)
//
// Setups are registered with the REGISTER_SETUP macro, which must be called
// once per setup in the user's setup.cpp file.
class SetupRegistry {
 public:
    using SetupFactory = std::function<std::unique_ptr<SetupBase>(
        Input&, Grid&, DataBlock&, Output&)>;

    /// Name reserved for the legacy Setup class (weak symbols)
    static constexpr const char* legacyName = "legacy";

    /// Register a setup under a given name (several setups may be registered)
    static void Register(const std::string &name, SetupFactory factory) {
        auto& instance = GetInstance();
        if(name.empty()) {
            IDEFIX_ERROR("Setup registration with an empty name is not allowed.");
        }
        if(name == legacyName) {
            std::stringstream msg;
            msg << "The setup name \"" << legacyName << "\" is reserved for the legacy "
                << "Setup class. Please choose another name." << std::endl;
            IDEFIX_ERROR(msg);
        }
        if(instance.factories.find(name) != instance.factories.end()) {
            std::stringstream msg;
            msg << "Duplicate registration of setup \"" << name << "\"." << std::endl
                << "Each setup must be registered under a unique name." << std::endl
                << "Registered setups are: " << GetStringOfNames() << std::endl;
            IDEFIX_ERROR(msg);
        }
        instance.factories.emplace(name, std::move(factory));
    }

    /// Create the setup selected in the input file ([Setup] name)
    static std::unique_ptr<SetupBase> Create(
        Input& input, Grid& grid, DataBlock& data, Output& output) {
        // Is a name requested in the input file?
        if(input.CheckEntry("Setup", "name") < 0) {
            auto& instance = GetInstance();
            if(!instance.factories.empty()) {
                idfx::cout << "Main: WARNING! No [Setup] name given in the input file, "
                           << "falling back to the legacy Setup." << std::endl
                           << "Main: Registered setups are: " << GetStringOfNames()
                           << std::endl;
            }
            return std::make_unique<DefaultSetup>(input, grid, data, output);
        }
        // A name is requested: use it
        std::string name = input.Get<std::string>("Setup", "name", 0);
        return Create(name, input, grid, data, output);
    }

    /// Create the setup with an explicit name ("legacy" for the legacy Setup)
    static std::unique_ptr<SetupBase> Create(const std::string &name,
        Input& input, Grid& grid, DataBlock& data, Output& output) {
        if(name == legacyName) {
            idfx::cout << "Main: using the legacy Setup." << std::endl;
            return std::make_unique<DefaultSetup>(input, grid, data, output);
        }
        auto& instance = GetInstance();
        auto it = instance.factories.find(name);
        if(it == instance.factories.end()) {
            std::stringstream msg;
            msg << "Setup \"" << name
                << "\" requested in the input file ([Setup] name) is not available." << std::endl
                << "Available setups are: " << GetStringOfNames() << std::endl
                << "Set [Setup] name to one of these, or to \"" << legacyName
                << "\" to use the legacy Setup class." << std::endl;
            IDEFIX_ERROR(msg);
        }
        idfx::cout << "Main: using setup \"" << name << "\"." << std::endl;
        return it->second(input, grid, data, output);
    }

    /// Names of the registered setups (sorted)
    static std::vector<std::string> GetRegisteredNames() {
        std::vector<std::string> names;
        for(const auto &entry : GetInstance().factories) {
            names.push_back(entry.first);
        }
        return names;
    }

 private:
    std::map<std::string, SetupFactory> factories;

    static SetupRegistry& GetInstance() {
        static SetupRegistry instance;  // Meyer's singleton
        return instance;
    }

    // Human-readable list of the registered setups (plus the legacy one)
    static std::string GetStringOfNames() {
        std::string list;
        for(const auto &name : GetRegisteredNames()) {
            if(!list.empty()) list += ", ";
            list += name;
        }
        if(!list.empty()) list += ", ";
        list += legacyName;
        return list;
    }

    SetupRegistry() = default;
    ~SetupRegistry() = default;
    SetupRegistry(const SetupRegistry&) = delete;
    SetupRegistry& operator=(const SetupRegistry&) = delete;
};

// Macros for one-line registration in user's setup.cpp
//
//  REGISTER_SETUP(MySetup)                  -> setup named "MySetup"
//  REGISTER_SETUP("MachReflection", MySetup) -> setup named "MachReflection"
//
// The name of the registrar object is built from the line number, so that
// several setups can be registered in the same file.
#define IDEFIX_SETUP_CONCAT_(a, b) a##b
#define IDEFIX_SETUP_CONCAT(a, b) IDEFIX_SETUP_CONCAT_(a, b)

// Pick the 1-argument or 2-argument form of REGISTER_SETUP
#define IDEFIX_SETUP_PICK_(_1, _2, NAME, ...) NAME

#define IDEFIX_SETUP_1_(SetupClass) \
    IDEFIX_SETUP_BODY_(#SetupClass, SetupClass, __LINE__)

#define IDEFIX_SETUP_2_(SetupName, SetupClass) \
    IDEFIX_SETUP_BODY_(SetupName, SetupClass, __LINE__)

#define IDEFIX_SETUP_BODY_(SetupName, SetupClass, Id) \
    namespace { \
        struct IDEFIX_SETUP_CONCAT(SetupRegistrar, Id) { \
            IDEFIX_SETUP_CONCAT(SetupRegistrar, Id)() { \
                SetupRegistry::Register(SetupName, \
                    [](Input& input, Grid& grid, DataBlock& data, Output& output) { \
                        return std::make_unique<SetupClass>( \
                            input, grid, data, output); \
                    }); \
            } \
        } IDEFIX_SETUP_CONCAT(setupRegistrar, Id); \
    }  // namespace

#define REGISTER_SETUP(...) \
    IDEFIX_SETUP_PICK_(__VA_ARGS__, IDEFIX_SETUP_2_, IDEFIX_SETUP_1_)(__VA_ARGS__)

#endif  // SETUP_REGISTRY_HPP_
