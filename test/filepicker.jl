@testset "filepicker.jl" begin

    function makeFixtureTree()
        tmp = mktempdir()
        mkpath(joinpath(tmp, "sub1"))
        mkpath(joinpath(tmp, "sub2"))
        mkpath(joinpath(tmp, ".hiddenDir"))
        write(joinpath(tmp, "a.txt"), "")
        write(joinpath(tmp, "b.zmx"), "")
        write(joinpath(tmp, ".hidden.txt"), "")
        write(joinpath(tmp, "sub1", "c.txt"), "")
        return tmp
    end

    function makeSortFixtureTree()
        tmp = mktempdir()
        mkpath(joinpath(tmp, "Zebra"))
        mkpath(joinpath(tmp, "apple"))
        write(joinpath(tmp, "b.png"), "")
        write(joinpath(tmp, "A.png"), "")
        write(joinpath(tmp, "c.txt"), "")
        write(joinpath(tmp, "D.md"), "")
        return tmp
    end

    @testset "listDirectory" begin
        tmp = makeFixtureTree()
        dirs, files = OpticTrace.listDirectory(tmp)
        @test dirs == [".hiddenDir", "sub1", "sub2"]
        @test files == [".hidden.txt", "a.txt", "b.zmx"]

        dirs2, files2 = OpticTrace.listDirectory(tmp; extensions = [".zmx"])
        @test dirs2 == [".hiddenDir", "sub1", "sub2"]
        @test files2 == ["b.zmx"]

        # only one level deep -- sub1's file isn't listed
        @test !("c.txt" in files)
        subDirs, subFiles = OpticTrace.listDirectory(joinpath(tmp, "sub1"))
        @test isempty(subDirs)
        @test subFiles == ["c.txt"]

        # includeHidden=false filters out dotfiles/dot-directories
        dirs3, files3 = OpticTrace.listDirectory(tmp; includeHidden = false)
        @test dirs3 == ["sub1", "sub2"]
        @test files3 == ["a.txt", "b.zmx"]

        @test_throws ArgumentError OpticTrace.listDirectory(tmp; sortBy = :bogus)
    end

    @testset "listDirectory sorting" begin
        tmp = makeSortFixtureTree()

        # default sortBy = :name is alphabetical and case-insensitive
        dirs, files = OpticTrace.listDirectory(tmp)
        @test dirs == ["apple", "Zebra"]
        @test files == ["A.png", "b.png", "c.txt", "D.md"]

        # sortBy = :type orders files by extension (case-insensitive), name as tie-breaker;
        # directories have no extension, so they still end up name-sorted among themselves
        dirsByType, filesByType = OpticTrace.listDirectory(tmp; sortBy = :type)
        @test dirsByType == ["apple", "Zebra"]
        @test filesByType == ["D.md", "A.png", "b.png", "c.txt"]
    end

    @testset "checkFilePickerSortMode" begin
        for sortBy in OpticTrace.FILE_PICKER_SORT_MODES
            @test OpticTrace.checkFilePickerSortMode(sortBy) === nothing
        end
        @test_throws ArgumentError OpticTrace.checkFilePickerSortMode(:bogus)
    end

    @testset "pathBreadcrumbs" begin
        crumbs = OpticTrace.pathBreadcrumbs("/Users/matt/data")
        @test crumbs == ["/" => "/", "Users" => "/Users", "matt" => "/Users/matt",
                          "data" => "/Users/matt/data"]

        @test OpticTrace.pathBreadcrumbs("/") == ["/" => "/"]

        # relative paths are normalized to absolute before splitting
        relCrumbs = OpticTrace.pathBreadcrumbs("relative/path")
        @test first(relCrumbs) == ("/" => "/")
        @test last(relCrumbs) == ("path" => abspath("relative/path"))
    end

    @testset "checkFilePickerMode" begin
        for mode in OpticTrace.FILE_PICKER_MODES
            @test OpticTrace.checkFilePickerMode(mode) === nothing
        end
        @test_throws ArgumentError OpticTrace.checkFilePickerMode(:bogus)
    end

    @testset "filePicker argument errors" begin
        tmp = makeFixtureTree()
        @test_throws ArgumentError OpticTrace.filePicker(tmp; mode = :bogus)
        @test_throws ArgumentError OpticTrace.filePicker(tmp; sortBy = :bogus)
        @test_throws ErrorException OpticTrace.filePicker(joinpath(tmp, "a.txt"))
    end

    @testset "Bonito UI smoke tests" begin
        import Bonito

        function renderToString(app)
            path = joinpath(mktempdir(), "out.html")
            Bonito.export_static(path, app)
            return read(path, String)
        end

        @testset "filePickerApp :file mode" begin
            tmp = makeFixtureTree()
            app = OpticTrace.filePickerApp(tmp; mode = :file)
            @test app isa Bonito.App
            html = renderToString(app)
            @test occursin("sub1", html)
            @test occursin("sub2", html)
            @test occursin("a.txt", html)
            @test occursin("b.zmx", html)
            @test occursin("Select", html)
            @test occursin("Cancel", html)
            @test occursin("ondblclick", html)
            @test occursin("Hidden", html)
            @test occursin("Name", html)
            @test occursin("Type", html)
            @test occursin("<select", html)
        end

        @testset "filePickerApp :directory mode" begin
            tmp = makeFixtureTree()
            app = OpticTrace.filePickerApp(tmp; mode = :directory)
            html = renderToString(app)
            @test occursin("sub1", html)
            # files are shown non-interactively for context in :directory mode
            @test occursin("a.txt", html)
            @test occursin("Select", html)
        end

        @testset "filePickerApp :multipleFiles mode" begin
            tmp = makeFixtureTree()
            app = OpticTrace.filePickerApp(tmp; mode = :multipleFiles)
            html = renderToString(app)
            @test occursin("a.txt", html)
            @test occursin("b.zmx", html)
            @test occursin("checkbox", html)
        end

        @testset "filePickerApp with extensions filter" begin
            tmp = makeFixtureTree()
            app = OpticTrace.filePickerApp(tmp; mode = :file, extensions = [".zmx"])
            html = renderToString(app)
            @test occursin("b.zmx", html)
            @test !occursin("a.txt", html)
        end

        @testset "filePickerApp hidden-file filtering" begin
            tmp = makeFixtureTree()

            # default (showHidden omitted) hides dotfiles/dot-directories
            appDefault = OpticTrace.filePickerApp(tmp; mode = :file)
            htmlDefault = renderToString(appDefault)
            @test !occursin(".hidden.txt", htmlDefault)
            @test !occursin(".hiddenDir", htmlDefault)
            @test occursin("a.txt", htmlDefault)

            # showHidden=true shows them
            appShown = OpticTrace.filePickerApp(tmp; mode = :file, showHidden = true)
            htmlShown = renderToString(appShown)
            @test occursin(".hidden.txt", htmlShown)
            @test occursin(".hiddenDir", htmlShown)
        end

        @testset "filePickerApp unreadable directory" begin
            if Sys.iswindows()
                @info "Skipping unreadable-directory test on Windows (chmod doesn't model POSIX permissions)"
            else
                tmp = mktempdir()
                blocked = joinpath(tmp, "blocked")
                mkpath(blocked)
                write(joinpath(blocked, "secret.txt"), "")
                chmod(blocked, 0o000)
                try
                    readable = try
                        readdir(blocked)
                        true
                    catch e
                        e isa Base.IOError ? false : rethrow()
                    end
                    if readable
                        @info "Skipping unreadable-directory test: readdir succeeded despite chmod 0o000 (likely running as root)"
                    else
                        # rendering rooted directly at the unreadable directory exercises the
                        # same code path as navigating into one via breadcrumbs/path field --
                        # the reactive callback runs on initial render too
                        app = OpticTrace.filePickerApp(blocked; mode = :file)
                        html = renderToString(app)
                        @test occursin("Cannot read this directory", html)
                        @test !occursin("secret.txt", html)
                    end
                finally
                    chmod(blocked, 0o755)
                end
            end
        end
    end

end
