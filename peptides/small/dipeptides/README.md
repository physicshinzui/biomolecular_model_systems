# Protonation-state variants

Amber/OpenMMの標準的なタンパク質残基テンプレートに対応する、ACE-X-NMEの
プロトン化状態バリアントです。

| ファイル | 状態 | 側鎖電荷 |
|---|---|---:|
| `ASH.pdb` | protonated Asp (`HD2` on `OD2`) | 0 |
| `GLH.pdb` | protonated Glu (`HE2` on `OE2`) | 0 |
| `HID.pdb` | neutral His, proton on `ND1` | 0 |
| `HIE.pdb` | neutral His, proton on `NE2` | 0 |
| `HIP.pdb` | doubly protonated His | +1 |
| `LYN.pdb` | neutral Lys | 0 |
| `CYM.pdb` | deprotonated Cys | -1 |

元のPDBファイルは変更していません。これらのファイルは次のコマンドで再生成できます。

```bash
python generate_protonation_states.py
```

OpenMMでは残基名を維持したまま読み込み、対応するAmber protein force fieldを
使用してください。`CYM`と`LYN`を含めOpenMMがACE-X-NME間のペプチド結合を確実に
認識できるよう、全共有結合を`CONECT`レコードにも明記しています。`CYX`は別の
CysとのS-S結合が必要なので、この単独残基セットには含めていません。
