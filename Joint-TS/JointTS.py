import torch
import torch.nn as nn
from layers.Transformer_EncDec import Encoder, EncoderLayer
from layers.SelfAttention_Family import FullAttention, AttentionLayer
from layers.Embed import PatchEmbedding
from utils.dlinear_utils import resolve_target_column


class _FlattenHead(nn.Module):
    def __init__(self, n_vars, nf, target_window, head_dropout=0):
        super().__init__()
        self.n_vars = n_vars
        self.flatten = nn.Flatten(start_dim=-2)
        self.linear = nn.Linear(nf, target_window)
        self.dropout = nn.Dropout(head_dropout)

    def forward(self, x):
        x = self.flatten(x)
        x = self.linear(x)
        x = self.dropout(x)
        return x


class _PatchBranch(nn.Module):
    def __init__(self, configs, n_vars, pred_len, patch_len=16, stride=8):
        super().__init__()
        self.seq_len = configs.seq_len
        self.pred_len = pred_len
        padding = stride

        self.patch_embedding = PatchEmbedding(
            configs.d_model, patch_len, stride, padding, configs.dropout)

        self.encoder = Encoder(
            [
                EncoderLayer(
                    AttentionLayer(
                        FullAttention(False, configs.factor, attention_dropout=configs.dropout,
                                      output_attention=False), configs.d_model, configs.n_heads),
                    configs.d_model,
                    configs.d_ff,
                    dropout=configs.dropout,
                    activation=configs.activation
                ) for _ in range(configs.e_layers)
            ],
            norm_layer=nn.Sequential(
                nn.LayerNorm(configs.d_model)
            )
        )

        self.head_nf = configs.d_model * int((configs.seq_len - patch_len) / stride + 2)
        self.head = _FlattenHead(n_vars, self.head_nf, pred_len, configs.dropout)

    def forward(self, x):
        means = x.mean(1, keepdim=True).detach()
        x = x - means
        stdev = torch.sqrt(torch.var(x, dim=1, keepdim=True, unbiased=False) + 1e-5)
        x /= stdev

        x = x.permute(0, 2, 1)
        enc_out, n_vars = self.patch_embedding(x)
        enc_out, _ = self.encoder(enc_out)
        enc_out = torch.reshape(enc_out, (-1, n_vars, enc_out.shape[-2], enc_out.shape[-1]))
        enc_out = enc_out.permute(0, 1, 3, 2)

        dec_out = self.head(enc_out)
        dec_out = dec_out.permute(0, 2, 1)

        dec_out = dec_out * (stdev[:, 0, :].unsqueeze(1).repeat(1, self.pred_len, 1))
        dec_out = dec_out + (means[:, 0, :].unsqueeze(1).repeat(1, self.pred_len, 1))
        return dec_out


class _GateEncoder(nn.Module):
    def __init__(self, d_model, n_heads, d_ff, e_layers, dropout):
        super().__init__()
        self.encoder = Encoder(
            [
                EncoderLayer(
                    AttentionLayer(
                        FullAttention(False, factor=1, attention_dropout=dropout, output_attention=False),
                        d_model, n_heads),
                    d_model,
                    d_ff,
                    dropout=dropout,
                    activation='gelu'
                ) for _ in range(max(1, e_layers))
            ],
            norm_layer=nn.Sequential(nn.LayerNorm(d_model))
        )
        self.proj = nn.Linear(d_model, 1)

    def forward(self, x):
                                                    
        enc_out, _ = self.encoder(x)
        gate = torch.sigmoid(self.proj(enc_out).squeeze(-1))
        return gate


class Model(nn.Module):
    def __init__(self, configs, patch_len=16, stride=8):
        super().__init__()
        self.pred_len = configs.pred_len
        target_idx, _ = resolve_target_column(configs.root_path, configs.data_path, configs.target)
        self.target_idx = max(0, min(target_idx, configs.enc_in - 1))
        self.surrogate_indices = [i for i in range(configs.enc_in) if i != self.target_idx]

        self.base_branch = _PatchBranch(configs, 1, configs.pred_len,
                                        patch_len=patch_len, stride=stride)
        sur_dim = max(1, len(self.surrogate_indices))
        self.sur_branch = _PatchBranch(configs, sur_dim, configs.pred_len,
                                       patch_len=patch_len, stride=stride)
        self.sur_proj = nn.Linear(sur_dim, 1)
                                               
        self.gate_encoder = _GateEncoder(d_model=1, n_heads=1, d_ff=16, e_layers=1, dropout=configs.dropout)

    def forward(self, x_enc, x_mark_enc, x_dec, x_mark_dec, mask=None):
        x_tar = x_enc[:, :, self.target_idx:self.target_idx + 1]
        base_pred = self.base_branch(x_tar)             

        if self.surrogate_indices:
            sur_input = x_enc[:, :, self.surrogate_indices]
        else:
            sur_input = torch.zeros(x_enc.size(0), x_enc.size(1), 1, device=x_enc.device, dtype=x_enc.dtype)
        sur_pred = self.sur_branch(sur_input)
        sur_residual = self.sur_proj(sur_pred).squeeze(-1)          

        gate = self.gate_encoder(base_pred)          
        residual = gate * sur_residual

        final = torch.zeros(base_pred.size(0), self.pred_len, x_enc.size(-1), device=base_pred.device,
                            dtype=base_pred.dtype)
        final[:, :, self.target_idx] = base_pred.squeeze(-1) + residual

        aux = {
            'base_target': base_pred.squeeze(-1),
            'residual': residual
        }
        return final, aux

    def regularization_loss(self):
        return None
