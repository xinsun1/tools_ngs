__author__ = 'sunxin'

'''
This script generate filemaker export list to zebra-label-printer pdf.
default label size = (35mm, 12mm)
dpi =
size = (350,120) pixel

For any questions, contact : sunxin931128@gmail.com
'''

import os, qrcode, sys

from PIL import Image, ImageDraw, ImageFont

out_size = (350,120)
label_size = (35,12)

qr = qrcode.QRCode(
    version=2,
    border=0
)

font_hel_bold = ImageFont.truetype('/Users/sunxin/Downloads/Arial-Unicode-Bold.ttf',
                                   size=20,
                                   index=0)

def generate_qr(id):

    global qr

    qr.clear()
    qr.add_data(str(id))
    qr_im = qr.make_image()

    return qr_im

def generate_label(id_list):

    global out_size,label_size, font_hel_bold

    out_im = Image.new("1", out_size, color=1)
    qr_im = generate_qr(id_list[0])
    out_im.paste(qr_im.resize((54,54)),(2,33))

    draw_out = ImageDraw.Draw(out_im)
    draw_out.text((60,4), str(id_list[1]).strip('"'), font=font_hel_bold)
    draw_out.text((60,32), str(id_list[2]).strip('"'), font=font_hel_bold)
    draw_out.text((60,60), str(id_list[3]).strip('"'), font=font_hel_bold)
    draw_out.text((60,88), str(id_list[4]).strip('"'), font=font_hel_bold)

    del draw_out

    return out_im

def print_label(file_h):
    fh = open(str(file_h), 'r')

    im_list = []
    while 1 :
        fl = fh.readline().strip().split(',')
        if len(fl) == 1 :
            break
            return 0

        im_list.append(generate_label(fl))

    im_list[0].convert('RGB').save('label.pdf',
                                   format='PDF',
                                   dpi=(254, 254),
                                   save_all=True,
                                   append_images=im_list[1:])


if __name__ == '__main__' :
    print('''usage : filemaker_qr_label COMMA_SEP_FILE

This script generate filemaker export list to zebra-label-printer pdf.
Default settings :
    sep = ','
    label size = (35mm, 12mm).
    dpi = (254, 254)
    size = (350,120) pixel
    output_file_name = label.pdf
    font_file = /Library/Fonts/Arial Unicode.ttf
    Require package : qrcode pillow
For any questions, contact : sunxin931128@gmail.com''')

    if len(sys.argv) == 1 :
        print()
        print('ERROR : no input file')
    else :
        print_label(sys.argv[1])



















