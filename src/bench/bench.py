"""
entry point for the spacer_bencher tool
"""
import sys


def main():
    """Main entry point that routes to appropriate command scripts"""
    if len(sys.argv) < 2:
        # No command provided, show help
        print("Spacer matching benchmark tool")
        print("Available commands:")
        print("  full_run        - Run the complete benchmarking pipeline")
        print("  simulate        - Generate simulated sequences")
        print("  generate_scripts - Generate tool execution scripts")
        print("  run_tools       - Execute tools on simulated data")
        print("  subsample       - Intelligently subsample real dataset")
        print("\nUse 'spacer_bencher <command> --help' for command-specific help")
        return
    
    command = sys.argv[1]
    
    # Route to appropriate command script
    command_scripts = {
        "full_run": "bench.commands.full_run",
        "simulate": "bench.commands.simulate", 
        "generate_scripts": "bench.commands.generate_scripts",
        "run_tools": "bench.commands.run_tools",
        "subsample": "bench.commands.subsample"
    }
    
    if command in command_scripts:
        # Remove the command from sys.argv and pass the rest to the subcommand
        sys.argv = sys.argv[1:]  # Remove the command name
        module_name = command_scripts[command]
        
        try:
            # Import and run the command module
            import importlib
            module = importlib.import_module(module_name)
            module.main()
        except ImportError as e:
            print(f"Error importing {module_name}: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"Error running {command}: {e}")
            sys.exit(1)
    else:
        print(f"Unknown command: {command}")
        print("Available commands: full_run, simulate, generate_scripts, run_tools, subsample")
        print("Use 'spacer_bencher <command> --help' for command-specific help")
        sys.exit(1)


if __name__ == "__main__":
    main()
