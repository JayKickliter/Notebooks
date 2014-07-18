using  DataFrames
import IPPDSP
import DSP
# import Base: isapprox


function polyize{T}( h::Vector{T}, interpolation )
    hLen         = length( h )
    tapsPerPhase = int( ceil( hLen/interpolation ))
    pfbSize = tapsPerPhase * interpolation
    # check that the vector is an integer multiple of interpolation
    if hLen != pfbSize 
        hExtended             = similar( h, pfbSize )
        hExtended[1:hLen]     = h
        hExtended[hLen+1:end] = 0
        h                     = hExtended
    end
    nFilters     = interpolation
    hLen         = length( h )
    tapsPerPhase = int( hLen/nFilters )
    pfb          = reshape( h, nFilters, tapsPerPhase )'
end

function upsample{T}( x::Vector{T}, interpolation )
    buffer = zeros(T, length(x) * interpolation )
    
    for n = 1:length(x)
        buffer[ n*interpolation - interpolation + 1 ] = x[n]
    end
    
    buffer
end

function interpolate{Th, Tx}( h::Vector{Th}, x::Vector{Tx}, interpolation )
    x       = upsample( x, interpolation )
    DSP.filt( h, x )
end

function polyinterpolate{T}( PFB::Array{T}, x::Vector{T} )
    (Ntaps, Nφ) = size( PFB )               # each column is a phase of the PFB, the rows hold the individual taps
    Nx          = length( x )               # number of input items    
    output      = similar( x, Nx * Nφ )     # Noutput = Nx * Nφ
    
    # until Xn == Ntaps, x[Xn-m+1] would reach out of bounds 
    #     this first loop limits Xn-m+1 to a minimum of 1
    #
    # for each input sample, for each φ (phase, or column of the PFB)
    for Xn = 1:Ntaps-1, φ = 1:Nφ                      
        output[Nφ*(Xn-1)+φ] = zero(T)
        # for each tap in phase[n]
        for Tn = 1:Xn
            @inbounds output[Nφ*(Xn-1)+φ] += PFB[Tn, φ] * x[Xn-Tn+1]
        end
    end
    
    # no longer in danger of stepping out of bounds
    #
    # for each input sample, φ (phase, or column of the PFB)
    for Xn = Ntaps:Nx, φ = 1:Nφ
        # for each tap in phase[n]
        output[Nφ*(Xn-1)+φ] = zero(T)
        for Tn = 1:Ntaps
            @inbounds output[Nφ*(Xn-1)+φ] += PFB[Tn, φ] * x[Xn-Tn+1]
        end
    end

    return output
end

polyinterpolate( h, x, interpolation ) = polyinterpolate( polyize( h, interpolation ), x )
polyinterpolate( h, x, interpolation ) = polyinterpolate( polyize( h, interpolation ), x )



function isapprox( x1::Vector, x2::Vector )
    Nx1 = length( x1 )
    Nx2 = length( x2 )    
    
    Nx1 == Nx2 || throw( DimensionMismatch("x1 & x2 must be the equal length vectors") )
    
    for i = 1:Nx1
        isapprox( x1[i], x2[i] ) || print( "Something when wrong at index ", i )
    end
end


function testspeed( Nsamples, Ntaps, interpolation )

    nativeTimes = zeros( length(Nsamples)*length(Ntaps) )
    ippTimes    = zeros( length(Nsamples)*length(Ntaps) )
    
    index      = 0

    for nt in Ntaps

    h = IPPDSP.lowpass( Float64, nt, 0.5/interpolation, true )

        for nx in Nsamples
            index += 1
            nx = int(nx)
            x  = rand( typeof(h[1]), nx )        
            
            gc()
            t = time()
            polyAns = polyinterpolate( h, x, interpolation )
            tPoly = time()-t
            nativeTimes[index] = tPoly
                    
            gc()
            t = time()
            ippAns = IPPDSP.filt( h, x, interpolation, 1 )
            tIPP = time()-t
            ippTimes[index] = tIPP
            
            areApprox( polyAns, ippAns )
        end
    end

    nativeDF = DataFrame( Size=numSamples, Time=nativeTimes, Implementation="Native" )
    ippDF    = DataFrame( Size=numSamples, Time=ippTimes, Implementation="IPP" )
    timeData = vcat( nativeDF, ippDF )     
    
end


function testaccuracy()
    interpolation = 4
    nx            = int(1e6) #100_000
    nt            = 46
    # round up so number of taps is an integer multiple of interpolation
    nt            = ceil(nt/interpolation)*interpolation
    
    x  = Float64[1:nx]
    h  = Float64[1:nt]
    # h  = IPPDSP.lowpass( Float64, nt, 0.5/interpolation, true )
        
    # gc()
    t = time()
    polyAns = polyinterpolate( h, x, interpolation )
    tPoly = time()-t

    # gc()
    t = time()
    ippAns = IPPDSP.filt( h, x, interpolation, 1 )
    tIPP = time()-t
    
    # display([ polyAns ippAns ])
    
    isapprox( polyAns, ippAns )
    
    @printf( "Poly Time: %0.3e\tIPP Time: %0.3e",  tPoly, tIPP )
end

times = testspeed(logspace(2,7, 2), 56, 4)
# testaccuracy()
