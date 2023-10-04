close all;
clear all;
clc;
%%project scaled version
%% Initial Design of the Filter (FIR), using Minimax Optimal
%specifications of the band-stop filter (notch filter)
wp1 = 0.15*pi; %passband 1 of the filter
wp2 = 0.35*pi; %passband 2 of the filter
ws1 = 0.19*pi; %stopband 1 of the filter
ws2 = 0.30*pi; %stopband 2 of the filter
M = 120; % order of the filter 
Stop_Atten = 40; %stopband attentuation in db
Pass_ripple = 0.05; %passband ripple in db
dp = (10^((Pass_ripple/2)/20)) - 1; %upper bound of passband ripple in non-db
ds = (10^(-Stop_Atten/20)); %stopband attentuation in non-db
% dp=0.0029 and ds=0.01 Therefore dp is more stringent
%setting up the minimax optimum filter
weight = [ds/dp, 1, ds/dp]; %put more weight on the passband
w = [0, wp1, ws1, ws2, wp2, pi]/pi; %setting the freq vector (w*pi) 
A = [1, 1, 0, 0, 1, 1]; %amplitude corresponding to the freq vector
Bstop1  = firpm(M, w, A, weight); %%creating the filter with park-mc


%% poles and zeros
[zero,pole,k] = tf2zpk(Bstop1, 1); %retuning zeros, poles, and gain constant

%reorder the zeros so that it is easy to group them
    zero_sort = sort(zero,'ComparisonMethod','real');

 
%% configuring the Biquads from the poles and zeros
    % complex conjugate coefficients result in real coeffcients
    % z = a +jb and z* = a - jb
    %each biquad will be in y[n] =  x[n] - (2a)x[n-1] + (a^2+b^2)x[n-2]
        biquad_new = zeros(M/2, 3);
        stage = 1;
            k_gain_stage = nthroot(k, max(length(biquad_new)));
    for zero_counter = 1:2:(M-1)
        %forming each biquads coeffcients
        biquad_new(stage, 1:3) = k_gain_stage*[1, (-2*real(zero_sort(zero_counter))), (real(zero_sort(zero_counter)))^2 + (imag(zero_sort(zero_counter)))^2];
        stage = stage + 1;
    end
    

    %combining two biquads together to form a quad (4 zeros/poles)
        quad = zeros(max(length(biquad_new))/2, 5);
        stage = 1;
       for stg_count = 1:1:(max(length(biquad_new))/2)
            quad(stage, :) = conv(biquad_new(stg_count, :), biquad_new(max(length(biquad_new))-stg_count+1,:)); 
            stage = stage + 1;
       end
        og_quad = quad;
        
      %combining two quads together (octoquad - 8 zeros/poles)
        octoquad = zeros(max(length(quad))/2, 9);
        stage = 1;
       for stg_count = 1:1:(max(length(quad))/2)
            octoquad(stage, :) = conv(quad(stg_count, :), quad(max(length(quad))-stg_count+1,:)); 
            stage = stage + 1;
       end


%% including the scaling factor
        
        
        %want to reorder the filters and have the most peaking at the end
        filt_gain = zeros(1, max(size(octoquad)));
        for stg = 1:1:max(size(octoquad))
            filt_gain(1, stg) = filternorm(octoquad(stg, :), 1, inf); %seeing ones with highest gain
        end
        
        for num1 = 1:1:(length(filt_gain)-1)
            if filt_gain(1, num1) > filt_gain(1, num1+1)
               filt_gain(1, [num1, num1+1]) = filt_gain(1, [num1+1, num1]);
               octoquad([num1+1, num1], :) = octoquad([num1, num1+1], :); %swap biquads (rows in matrix)
               loop = true;
               num2 = num1+1;
                if(num1 < length(filt_gain)-1)
                   while loop == true
                        if filt_gain(1, num2) > filt_gain(1, num2+1)
                            filt_gain(1, [num2, num2+1]) = filt_gain(1, [num2+1, num2]);
                            octoquad([num1+1, num1], :) = octoquad([num1, num1+1], :); %swap biquads (rows in matrix)
                            num2 = num2+1;
                            if num2 == length(filt_gain)
                                loop = false;
                            end                           
                        else
                            loop = false;
                        end
                   end
                end 
            end
        end
        for num1 = 2:1:(length(filt_gain))
            if filt_gain(1, num1) < filt_gain(1, num1-1)
               filt_gain(1, [num1, num1-1]) = filt_gain(1, [num1-1, num1]);
               octoquad([num1-1, num1], :) = octoquad([num1, num1-1], :); %swap biquads (rows in matrix)
               loop = true;
               num2 = num1-1;
               if(num2 > 1)
                   while loop == true
                        if filt_gain(1, num2) < filt_gain(1, num2-1)
                            filt_gain(1, [num2, num2-1]) = filt_gain(1, [num2-1, num2]);
                            octoquad([num2-1, num2], :) = octoquad([num2, num2-1], :); %swap biquads (rows in matrix)
                            num2 = num2-1;
                            if num2 == 1
                                loop = false;
                            end 
                        else
                            loop = false;
                        end
                   end
               end 
            end
        end        
 
 
   %F filter relates input to intermediate node (v[n])
   %{
           We will not find F for the last stage because that wont have any
           scaling as the output will undo all the scaling done prev
   %}
    F(1, :) = octoquad(1, :); %setting initial F filter
    for stg = 2:1:max(size(octoquad)-1) %finding each of the F filters
        F_stage_calc = conv(octoquad(stg, :), F((stg-1), :));
        F(stg, 1:length(F_stage_calc)) = F_stage_calc; 
    end
    %{
    ---R is the bound set for the F.T.(input)
    ---input in this case is a sinusoid where F.T. is two delta functionns
    with with a magnitude of 1/2 of the sinusoid amplitude
    
    ---our input is said to be a -6db below unity sinusoid power wise
    ---power of a sinusoid = (A^2)/2
    ---power of unity sinusoid = 10log(1/2)
    ---10log((A^2)/2) - 10log(1/2) = -6dB
    ---we find A = 0.5012
    %}
    R = 0.5012; %one norm of input for sinusoid when A = 0.5012
    
    %alpha exists when you need to scale initially for the first feedback
    % usually alpha = (1/R)/filternorm(F(1), 1, 2); for IIR filters 
    alpha = 1; %fir filters have alpha set as 1 since no feedback
    lambda_total = alpha; %setting initial lambda
    % only goes up to second last stage because last lambda undoes scaling
    for stg = 1:1:max(size(octoquad)-1)
        lambda(stg) = (1/R)/(lambda_total*filternorm(F(stg, :), 1, inf));
        %lambda total = the multplcation of all the scaling factors so far
        lambda_total = lambda_total*lambda(stg);
    end
        lambda(max(size(octoquad))) = 1/lambda_total; %undoing all scaling for last stage
    
    %applying lambda to each of the biquads
    for stg = 1:1:max(size(octoquad))
        octoquad(stg,:) = lambda(stg) * octoquad(stg,:); 
    end
    
    %quantizing coefficeints of the each stage
    for quad_cnt = 1:1:max(size(octoquad))
        octoquad_hold = octoquad(quad_cnt,:);
        octoquad_quant = quantizing(octoquad_hold, 3, 14);
        octoquad(quad_cnt,:) = octoquad_quant;
    end

    %checking the octoquads come together properly
            octoquad1 = octoquad(1,:);
            %convolve all quad stages together
        for stg_count = 1:1:(max(length(octoquad))-1)
            octoquad2 = octoquad(stg_count+1,:);
            full_octoquad = conv(octoquad1, octoquad2);
            octoquad1 = full_octoquad;
        end

%% using stages from the previous section to test input/output (time domain)
    %below is different frequencies to be tested
    Fs = 10^6; %(1MHz)
    omega1 = 95000*2*pi; %Hz * 2pi (rad/sec) of first sinusoid
    omega2 = 50000*2*pi; %Hz * 2pi (rad/sec) of second sinusoid
    omega3 = 200000*2*pi; %Hz * 2pi (rad/sec) of third sinusoid
    A = 0.5012; %ampltiude of the input signal
    time = 0:1/(Fs*100):4*10^-4;
    Ts = 0:1/Fs:max(time);

   %Quantizing adding dither to input signal and then quantizing it 
    %{
     rand function is uniformly distributed while randn is normally
    distributed
    %}
    %we want a range from [-delta/2 to +delta/2]
    %rand generates a number from 0 to 1 uniformly
    %{
    also note that diither_b must equal the bits used to quantized the
    signal in order to ensure harmonics are removed
    %}
    dither_b = 15;
    delta = 2^-dither_b;
    %want to shift the rand function from [0 1] to [-1 1]
    dither = (delta/2)*((rand(1, length(Ts))-0.5)*2);

 %first input
    input_1 = A*sin(omega2*Ts); %setting input
    dither_input_1 = input_1 + dither; %adding dither to quantize signal
    quant_dither_input_1 = quantizing(dither_input_1,1,15); %1 bit for sign and 11 frac bits
    test_input_1 = quant_dither_input_1;

    dec_bits = 2;   % decimal + sign bits for coefficients of filter stages
    frac_bits = 23; % fraction bits for coefficients of filter stages
    %calling function that implements input through filter
    matrix_output_1 = filtering(test_input_1, octoquad, dec_bits, frac_bits);
    dim_matrix_output_1 = size(matrix_output_1); %dim of matrix_output for following line
    %want the final row of matrix_output
    test_output_1 = matrix_output_1(dim_matrix_output_1(1,1) , :);
    

 %second input
    input_2 = A*sin(omega3*Ts); %setting input
    dither_input_2 = input_2 + dither; %adding dither to quantize signal
    quant_dither_input_2 = quantizing(dither_input_2,1,15); %1 bit for sign and 11 frac bits
    test_input_2 = quant_dither_input_2;

    %calling function that implements input through filter
    matrix_output_2 = filtering(test_input_2, octoquad, dec_bits, frac_bits);
    dim_matrix_output_2 = size(matrix_output_2); %dim of matrix_output for following line
    %want the final row of matrix_output
    test_output_2 = matrix_output_2(dim_matrix_output_2(1,1) , :);

 %third input
    input_3 = A*sin(omega1*Ts); %setting input
    dither_input_3 = input_3 + dither; %adding dither to quantize signal
    quant_dither_input_3 = quantizing(dither_input_3,1,15); %1 bit for sign and 11 frac bits
    test_input_3 = quant_dither_input_3;
    %calling function that implements input through filter
    matrix_output_3 = filtering(test_input_3, octoquad, dec_bits, frac_bits);
    dim_matrix_output_3 = size(matrix_output_3); %dim of matrix_output for following line
    %want the final row of matrix_output
    test_output_3 = matrix_output_3(dim_matrix_output_3(1,1) , :);



%% Plots
    %plotting impulse both impulse response (before and after quantization)
    figure(1)
        subplot(3,1,1)
        stem(Bstop1);
        grid on;
        title('Impulse Reponse Non-Quantized');
        xlabel('Samples');
        ylabel('Amplitude');
    
        subplot(3,1,2);
        stem(full_octoquad);
        grid on;
        title('Impulse Reponse with Quantized Coefficients in the Stages');
        xlabel('Samples');
        ylabel('Amplitude');

        subplot(3,1,3);
        stem(Bstop1, 'x');
        hold on;
        stem(full_octoquad);
        title('Non-Quantized vs Quantized Impulse Response');
        legend('Non-Quantized', 'Quantized');
        xlabel('Samples');
        ylabel('Amplitude');
        grid on;
    
    %Plotting the freq-response of our non-quant and quant filter
      figure(2)
        %non quantized Bandstop filter
        [H_bandstop1_quant, freq_axis_quant] = freqz(Bstop1, 1, 1024, 'Whole');
        subplot(3, 1, 1);
        plot(freq_axis_quant/pi, 20*log10(abs(H_bandstop1_quant)));
        grid on;
        title('Non Quantized Bandstop Freq. Response');
        xlabel('Normalized Frequency (x\pi rad/sample)')
        ylabel('Magnitude (dB)')
    
        %quantized Bandstop Filter
        full_octoquad_zeropad = horzcat(full_octoquad, zeros(1,500));
        N_DFT_full_octoquad = length(full_octoquad_zeropad);
        freq_res_full_octoquad = 2*pi/N_DFT_full_octoquad;
        freq_bin_full_octoquad = (0:1:N_DFT_full_octoquad-1)*freq_res_full_octoquad;
        dft_full_octoquad = fft(full_octoquad_zeropad);
        subplot(3,1,2)
        plot(freq_bin_full_octoquad/pi, 20*log10(abs(dft_full_octoquad)));
        grid on;
        title('Quantized Bandstop Freq. Response');
        xlabel('Normalized Frequency (x\pi rad/sample)')
        ylabel('Magnitude (dB)')
        
        subplot(3,1,3)
        plot(freq_axis_quant/pi, 20*log10(abs(H_bandstop1_quant)), '-x');
        hold on;
        plot(freq_bin_full_octoquad/pi, 20*log10(abs(dft_full_octoquad)), '-o');
        grid on;
        title('Non-Quantized vs Quantized Freq. Response');
        xlabel('Normalized Frequency (x\pi rad/sample)')
        ylabel('Magnitude (dB)')
        legend('Non-Quantized', 'Quantized');
        
    
   %%plotting input and output tests
    %plotting the first input and outpuot
    figure(3)
        subplot(3,1,1);
        plot(Ts, test_input_1);
        grid on;
        title('Input Signal 1 of 50Khz')
        xlabel('time (sec)');
        ylabel('Amplitude (Voltage)')

        subplot(3,1,2)
        plot(Ts, test_output_1);
        grid on;
        title('Output Signal 1 of 50Khz');
        xlabel('time (sec)');
        ylabel('Amplitude (Voltage)')


        subplot(3,1,3)
        plot(Ts, test_input_1, '-x');
        grid on;
        hold on;
        plot(Ts, test_output_1, '-o');
        title('Input vs Output Signal 1 of 50Khz');
        xlabel('time (sec)');
        ylabel('Amplitude (Voltage)')
        legend('Input 1', 'Output 1');
    
    
    %%plotting the Second input and output
    figure(4)
        subplot(3,1,1);
        plot(Ts, test_input_2);
        grid on;
        title('Input Signal 2 of 200Khz')
        xlabel('time (sec)');
        ylabel('Amplitude (Voltage)')

        subplot(3,1,2)
        plot(Ts, test_output_2);
        grid on;
        title('Output Signal 2 of 200Khz');
        xlabel('time (sec)');
        ylabel('Amplitude (Voltage)')

        subplot(3,1,3)
        plot(Ts, test_input_2, '-x');
        grid on;
        hold on;
        plot(Ts, test_output_2, '-o');
        title('Input vs Output Signal of 200Khz');
        xlabel('time (sec)');
        ylabel('Amplitude (Voltage)')
        legend('Input 2', 'Output 2');


    %%plotting the Third input and output
    figure(5)
        subplot(3,1,1);
        plot(Ts, test_input_3);
        grid on;
        title('Input Signal 3 of 95Khz')
        xlabel('time (sec)');
        ylabel('Amplitude (Voltage)')
    
        subplot(3,1,2)
        plot(Ts, test_output_3);
        grid on;
        title('Output Signal 3 of 95Khz');
        xlabel('time (sec)');
        ylabel('Amplitude (Voltage)')
    
        subplot(3,1,3)
        plot(Ts, test_input_3, '-x');
        grid on;
        hold on;
        plot(Ts, test_output_3, '-o');
        title('Input vs Output Signal 3 of 95Khz');
        xlabel('time (sec)');
        ylabel('Amplitude (Voltage)')
        legend('Input 3', 'Output 3');

    
    %calculating and plotting SNR     
    figure(6)
    snr_function_of_output_1 = snr(test_output_1)
    snr(test_output_1);

    figure(7)
    snr_function_of_output_2 = snr(test_output_2)
    snr(test_output_2);

%% FVtool Plots
    %z-plane zeros and poles
    zeropole_plot = fvtool(Bstop1,1, full_octoquad, 1, "polezero");
    legend(zeropole_plot, "Non-Quantized", "Quantized");

    %Group Delay
    grp_delay = fvtool(Bstop1, 1, full_octoquad, 1, 'grpdelay');
    legend(grp_delay, 'Non-Quantized', 'Quantized')



%% Estimated SNR
dim_stages_n = size(octoquad); %find amound of stages and coefficients
num_stages = dim_stages_n(1,1);
num_coeff = dim_stages_n(1,2);

b = 23; %fraction bits
variance = ((2^-b)^2)/12;
noise = zeros(1,num_coeff);
for coeff_num = 1:1:num_stages
    noise(coeff_num) = num_coeff*variance;
end

% S(out) = S(in)*|H(w)|^2
%for total noise we can just integrate S(out) and divide over the period (2pi)
%causes |H(w)|^2 to become the the square of the sec norm (no sqrt)
%H(w) is the filter the noises see (e.g. first noise sees every stage except first


for g_num = 1:1:num_stages-1
    G_calc = octoquad(g_num+1, :);
    for stg_cnt = (g_num+1):1:num_stages-1
        G_calc = conv(octoquad(stg_cnt+1,:), G_calc(:));
    end 
    G(g_num, 1:length(G_calc)) = G_calc;
    clear G_calc;
end    

    tot_noise_out = 0;
    for stg_num = 1:1:num_stages-1
        tot_noise_out =  noise(stg_num)*(filternorm(G(stg_num,:),1,2))^2 + tot_noise_out;
    end
        tot_noise_out = tot_noise_out + noise(num_stages);

    
    signal_noise_ratio = 10*log10(((R^2)/2) / tot_noise_out);
    fprintf('My own Calculated SNR done in Matlab, assuming white noise, uniform Distribution and independent of each other: %4.2f dB\n',  signal_noise_ratio)




%% Rounding Quantization function using Quant
function quant_n = quantizing(pre_quant, dec_bits, frac_bits)
    quant_n = zeros(1, length(pre_quant));
    b_frac = frac_bits; %fraction bit length
    b_dec = dec_bits;
    max_dec = 0; %intial initialization of maximum bit size
    delta = 2^-b_frac; %step difference between quant levels

    for max_dec_calc = 0:1:(b_dec-1) %highest number in q should dec bits + fraction bits
        max_dec = 2^max_dec_calc + max_dec;
    end

    if b_dec > 1 && b_frac > 1
        max_num = max_dec + 1 - delta;
    else
        max_num = 1 - delta;
    end

    
    quant_n(1,:) = quant(pre_quant, delta);
    for num_pos = 1:1:length(quant_n)   %if number is greater than range then go to max number
        if abs(quant_n(num_pos)) >= max_num %if reaches highest number then round to it
            if quant_n(num_pos) > 0    %goes above the bounds
            quant_n(num_pos) = max_num;
            else                        %goes below the bounds
                quant_n(num_pos) = -max_num;
            end
        end
    end
end


%% function for going through stages, utilizes quantizing function from above aswell
function output_n = filtering(input_n, stages_n, dec_bits, frac_bits)
    dim_stages_n = size(stages_n); %find amound of stages and coefficients
    num_stages = dim_stages_n(1,1);
    num_coeff = dim_stages_n(1,2);
    coeff_n = zeros(size(stages_n));

    % getting all the coefficeints set in the right vector
    for stage_cnt = 1:1:num_stages
            coeff_n(stage_cnt, :) = stages_n(stage_cnt, :);
    end
    
    test_input_n = horzcat(zeros(1,num_coeff), input_n); %want to make sure there are 0's for time<0
    %going to have a matrix that has every output of each stage including the final output 
    matrix_input_n = zeros(num_stages+1, length(test_input_n)); 
    matrix_input_n(1,:) = test_input_n;

    for input_cnt = num_coeff:1:(length(test_input_n))
        for stage_cnt = 1:1:num_stages
            result_stage = 0;
            for coeff_cnt = 1:1:num_coeff %run through each tap for the stage
                % e.g. c0*x[n] + C1*x[n-1] + C2*x[n-2]      %round aswell
                result_stage = quantizing(matrix_input_n(stage_cnt, input_cnt-coeff_cnt+1)*coeff_n(stage_cnt, coeff_cnt), dec_bits, frac_bits) + result_stage;
                result_stage = quantizing(result_stage, dec_bits, frac_bits);
            end

        matrix_input_n(stage_cnt+1, input_cnt) = result_stage;
        end 
    end 
    output_n = matrix_input_n(:, (num_coeff+1):length(test_input_n));
end





