---
layout: gallery
permalink: /design/
title: design
category: design
nav: true
nav_order: 3
pagination:
    enabled: false
---


<div class="post">

{% assign blog_name_size_design = site.blog_name_design | size %}
{% assign blog_description_size_design = site.blog_description_design | size %}

{% if blog_name_size_design > 0 or blog_description_size_design > 0 %}

  <div class="header-bar">
    <h1>{{ site.blog_name_design }}</h1>
    <h2>{{ site.blog_description_design }}</h2>
  </div>
  {% endif %}

</div>


<div class="row row-cols-1 row-cols-sm-2 row-cols-md-3 g-4">
  {% assign design_images = site.static_files | where_exp: "file", "file.path contains '/assets/img/design/'" %}
  {% for image in design_images %}
  <div class="col">
    <div class="card h-100 shadow-sm">
      <a href="{{ image.path | relative_url }}" target="_blank">
        <img src="{{ image.path | relative_url }}" class="card-img-top" alt="design image" style="object-fit: cover; height: 250px;">
      </a>
    </div>
  </div>
  {% endfor %}
</div>

