function waveNumber(folder, nemohInput)
    freq = zeros(size(nemohInput.wavefreq));
    if nemohInput.wavetype==1
        freq=nemohInput.wavefreq;
    elseif nemohInput.wavetype==2
        freq=nemohInput.wavefreq/(2*pi);
    else 
        freq=(2*pi) ./ nemohInput.wavefreq;
    end
    k_numbers = WNumber(freq, nemohInput.depth);
    fid = fopen(fullfile(folder,'WaveNumber.dat'),'w');
    for ii = 1:length(k_numbers)
        fprintf(fid,'%g\t %g\t %g\n', freq(ii), k_numbers(ii), (2*pi)/freq(ii));
    end
    fclose(fid);
end

function k_numbers = WNumber(omega, h)
    g = 9.81; 
    omega = omega(:); 
    k_numbers = zeros(size(omega));
    for i = 1:length(omega)
        w = omega(i);
        k0 = w^2 / g; 
        k_numbers(i) = fzero(@(k) w^2 - g*k*tanh(k*h), k0);
    end
end
