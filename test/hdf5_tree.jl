using HDF5

@testset "HDF5 inspection and legacy I/O" begin
    mktempdir() do dir
        path = joinpath(dir, "workspace.nxs")
        h5open(path, "w") do file
            file["data/axes"] = [1, 2, 3]
            file["data/scalar"] = 42
            file["data/large"] = zeros(10_000)
            file["data/title"] = repeat("α", 2_000)
            attributes(file["data/axes"])["long_name"] = "[H,0,0]"
            for i in 1:201
                file["runs/experiment$i/instrument_xml"] = repeat("xml", 1_000)
            end
        end
        io = IOBuffer()
        h5tree(path; keypath="/data", attributes=true, values=true, maxstring=12, io)
        output = String(take!(io))
        @test occursin("@long_name = \"[H,0,0]\"", output)
        @test occursin("= [1, 2, 3]", output)
        @test occursin("= 42", output)
        @test occursin(repeat("α", 12) * "…", output)
        @test occursin(r"large +\[Float64: 10000\]\n", output)
        @test ncodeunits(output) < 1_000

        h5tree(path; keypath="/runs", values=true, maxchildren=2, maxstring=16, io)
        output = String(take!(io))
        @test occursin("199 more children", output)
        @test ncodeunits(output) < 1_000

        h5tree(path; keypath="data/axes", values=true, io)
        @test occursin("= [1, 2, 3]", String(take!(io)))
        h5tree(path; maxdepth=0, io)
        @test String(take!(io)) == path * "\n"
        @test_throws ArgumentError h5tree(path; maxchildren=-1, io)

        legacy = joinpath(dir, "legacy")
        result = LSWTResults(PHysicalTDA.Sunny.Intensities{Float64,PHysicalTDA.Sunny.QPath,2}[],
                             PHysicalTDA.Sunny.QPath[], [1.0]; meta=Dict{String,Any}("tag" => "legacy"))
        saveHDF5("fixture", "temperature", "5K", String[], result; fname=legacy)
        saveHDF5("fixture", "temperature", "10K", String[], result; fname=legacy)
        @test loadHDF5(legacy; keypath="crystals/fixture/LSWT/temperature/10K/meta/tag") == "legacy"
        @test loadHDF5(legacy; keypath="crystals/fixture/LSWT/temperature/5K/meta/tag") == "legacy"
        h5tree(legacy; io)
        @test startswith(String(take!(io)), legacy * ".h5\n")
    end
end
