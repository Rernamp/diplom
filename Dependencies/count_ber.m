function [ BER ] = count_ber( message_RX, message_TX)

   N_bits = length(message_TX);
   N_errors = sum(xor(message_RX, message_TX));
   
   BER = N_errors/N_bits; 
end

