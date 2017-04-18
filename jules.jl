# global variables

Fs = 5000; # sampling frequency
T = 1/Fs;
# functions

function gentime(trial, sampling_freq=Fs)
    T = 1/sampling_freq
    L = size(trial, 1)
    return t = collect(0:L-1)*T
end

function genrange(start, finish)
    return convert(Integer, floor(start)):convert(Integer, floor(finish))
end

function genrangeseconds(start_sec, finish_sec, sampling_freq=Fs)
    return genrange(start_sec*sampling_freq, finish_sec*sampling_freq)
end

function fftplot(X; sampling_freq=Fs)
    T = 1/sampling_freq
    L = size(X, 1)
    t = collect(0:L-1)*T

    Y = fft(X)
    P2 = abs(Y/L)
    P1 = P2[1:div(L,2)+1]
    P1[2:end - 1] = 2*P1[2:end-1]
    
    freq = sampling_freq*collect(0:(L/2))/L
    
    subplot(2,1,1)
    plot(t,X)
    title("Time domain")
    xlabel("time (s)")
    ylabel("Amplitude")
    subplots_adjust(hspace=.5)
    subplot(2,1,2)
    plot(freq,P1)
    title("Single-Sided Amplitude Spectrum")
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    return (P1,freq)
end

function subdc(trial)
    trial - mean(trial)
end

function discreteder(trial, sampling_freq=Fs)
    T_s = 1/Fs
    b = zeros(size(trial))

    for i in 1:length(trial)-1
        b[i] = (trial[i+1] - trial[i])/T_s
    end

    return b
end

function getstats(trial_column)
    stddev = std(trial_column)
    mn = mean(trial_column)
    println("min: ", minimum(trial_column))
    println("max: ", maximum(trial_column))
    println("average: ", mn)
    println("standard dev: ", stddev)
    ran = -3.5*stddev+mn : 0.01 : 3.5*stddev+mn
    #plot(ran, gaussian(ran, mn, stddev))
    #hold
    #plot(stddev+mn,0,"o")
    #plot(-stddev+mn,0,"o")
end

function gaussian(x, avg, std_dev)
    return 1/(std_dev * sqrt(2*pi))*exp(-.5*((x-avg)/std_dev).^2)
end

function smooth(arr, num_point=5)
    if num_point % 2 == 0
        println("num_point must be odd")
        return 
    end
        
    b = zeros(size(arr))
    span = div(num_point, 2)
    start, finish = 1, length(arr)
    
    for i in 1:length(arr)
        n = 0
        
        inspan = i-span < start ? i-start : span
        inspan = i+span > finish ? finish-i : inspan
        
        #println("start ", start, " finish ", finish, " i ", i, " inspan ", inspan)
        for j in max(start, i-inspan):min(finish, i+inspan)
            #println("j", j, " arr ", arr[j])
            b[i] += arr[j]
            n+=1
        end
        
        b[i]/=n
    end
    return b
end

#LOW-PASS WINDOWED-SINC FILTER

function windowsinc(input_arr, fc, kernel_length::Integer)
    if kernel_length % 2 == 0
        println("num_point must be odd")
        #only applies to 1 indexed languages
        return 
    end
    #Set the cutoff frequency (between 0 and 0.5 relative to sampling freq)
    #kernel_length = 4/BW

    input = vcat(zeros(kernel_length), input_arr)
    kernel = zeros(kernel_length) #holds the filter kernel
    output = zeros(length(input))

    #Calculate the low-pass filter kernel

    kernel_length_half = div(kernel_length+1, 2)
    for i in 1:kernel_length
            if (i-kernel_length_half) == 0 
                    kernel[i] = 2*pi*fc
            else
                    kernel[i] = sin(2*pi*fc * (i-kernel_length_half)) / (i-kernel_length_half)
                    kernel[i] = kernel[i] * (0.54 - 0.46*cos(2*pi*i/kernel_length) )
            end
    end


    #Normalize the low-pass filter kernel for unity gain at DC
    kernel_sum = sum(kernel)
    kernel = kernel ./ kernel_sum

    #Convolve the input signal & filter kernel
    for j in kernel_length+1:length(input) 
            for i in 1:kernel_length
                    output[j] += input[j-i] * kernel[i]
            end
    end
    return output[kernel_length+1:end]
end

function acceltogs(x)
    return x*2.0/32768
end

function windowsinchp(input_arr, fc, kernel_length::Integer)
    if kernel_length % 2 == 0
        println("num_point must be odd")
        #only applies to 1 indexed languages
        return 
    end
    #Set the cutoff frequency (between 0 and 0.5 relative to sampling freq)
    #kernel_length = 4/BW

    input = vcat(zeros(kernel_length), input_arr)
    kernel = zeros(kernel_length) #holds the filter kernel
    output = zeros(length(input))

    #Calculate the low-pass filter kernel

    kernel_length_half = div(kernel_length+1, 2)
    for i in 1:kernel_length
            if (i-kernel_length_half) == 0 
                    kernel[i] = 2*pi*fc
            else
                    kernel[i] = sin(2*pi*fc * (i-kernel_length_half)) / (i-kernel_length_half)
                    kernel[i] = kernel[i] * (0.54 - 0.46*cos(2*pi*i/kernel_length) )
            end
    end


    #Normalize the low-pass filter kernel for unity gain at DC
    kernel_sum = sum(kernel)
    kernel = kernel./kernel_sum

    kernel = kernel .* -1
    kernel[kernel_length_half] += 1
    
    #Convolve the input signal & filter kernel
    for j in kernel_length+1:length(input) 
            for i in 1:kernel_length
                    output[j] += input[j-i] * kernel[i]
            end
    end
    return output[kernel_length+1:end]
end 

function findpeaks(x, tolerance=0.0001)
    peaks = []
    counter::Integer = 1
    for i in 2:length(x)-1
        if abs(x[i] - x[i+1]) > tolerance && abs(x[i] - x[i-1]) > tolerance
            append!(peaks, [(x[i], i)])
            counter += 1
        end
    end
    return peaks
end

function decimate(x, num::Integer)
    return x[1:num:length(x)]
end

function genfreqindex(finish, sampling_freq = Fs)
    return 1:finish/sampling_freq
end

function getenergy(signal)
    return sum(signal.^2)
end

function nernst(extra, intra, z, T=298)
    R = 8.314472
    F = 9.64853399 * 10.0^4
    E = (R*T)/(z*F)log(extra/intra)
end

function idft(X)
    a = conj(X)
    b = fft(a)
    return real(b)
end
